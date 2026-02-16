"""
Feigenbaum bifurcation analysis for the resonance collapse dynamics.

The formalism predicts that the Feigenbaum constant δ ≈ 4.669 appears in the
collapse dynamics. We test this by varying the coupling strength γ and detecting
period-doubling bifurcations in the steady-state oscillation of the resonance
expectation value ⟨R̂⟩.

Additionally, we test the spectral determinant prediction:
    det(R̂ - λI) = Π_{p≤N} (1 - λ/p(p-1))
which stabilizes at λ ≈ 0.002.
"""

from typing import List, Tuple, Dict
import numpy as np
from .prime_utils import sieve_primes
from .collapse_hamiltonian import (
    build_full_hamiltonian,
    build_resonance_operator,
    run_collapse_simulation,
    compute_entropy,
    compute_resonance_expectation,
)


def spectral_determinant(primes: List[int], lam: float) -> float:
    """
    Compute the spectral determinant from the formalism:

    det(R̂ - λI) = Π_{p≤N} (1 - λ / (p(p-1)))

    This product converges as N → ∞ and its zero defines the
    symbolic eigenvalue λ*.

    Returns:
        Product value (approaches 0 at λ*)
    """
    product = 1.0
    for p in primes:
        if p < 2:
            continue
        factor = 1.0 - lam / (p * (p - 1))
        product *= factor
    return product


def find_spectral_eigenvalue(
    primes: List[int],
    lam_min: float = 0.0001,
    lam_max: float = 1.0,
    n_points: int = 10000
) -> Tuple[float, List[Tuple[float, float]]]:
    """
    Find the smallest positive eigenvalue where the spectral determinant
    approaches its minimum (zeros correspond to resonance eigenvalues).

    Returns:
        (eigenvalue, [(lambda, det_value), ...])
    """
    lambdas = np.linspace(lam_min, lam_max, n_points)
    values = [(float(l), spectral_determinant(primes, float(l))) for l in lambdas]

    # Find the lambda where det is closest to zero (smallest positive root)
    min_abs = float('inf')
    best_lambda = 0.0
    for l, v in values:
        if abs(v) < min_abs and v > 0:
            # We want the value just before it crosses zero
            min_abs = abs(v)
            best_lambda = l

    return best_lambda, values


def compute_eigenvalue_ratios(
    primes: List[int],
    gamma: float = 1.0,
    use_legendre: bool = True
) -> Dict[str, object]:
    """
    Compute eigenvalue ratios of the full Hamiltonian and compare to
    physical constant ratios predicted by the twist framework.

    Key predictions tested:
    - 17 (trefoil complexity) in eigenvalue ratios
    - 137 (fine structure) in eigenvalue spacings
    - 108 (twist unit) in eigenvalue sums/products

    Returns:
        Dict with eigenvalues, ratios, and matching scores.
    """
    H = build_full_hamiltonian(primes, gamma, use_legendre)
    eigenvalues = np.linalg.eigvals(H)

    # Sort by real part
    sorted_eigs = sorted(eigenvalues, key=lambda e: e.real)
    real_parts = np.array([e.real for e in sorted_eigs])

    # Compute consecutive ratios
    ratios = []
    for i in range(len(real_parts) - 1):
        if abs(real_parts[i]) > 1e-10:
            ratios.append(abs(real_parts[i + 1] / real_parts[i]))

    # Compute spacings
    spacings = np.diff(real_parts)

    # Look for key numbers in ratios and spacings
    target_numbers = {
        17: "trefoil complexity (m_p/m_e factor)",
        29: "boundary prime (α correction)",
        108: "twist unit (2²×3³)",
        137: "fine structure constant",
    }

    matches = {}
    for target, meaning in target_numbers.items():
        # Check if any ratio is close to the target
        best_ratio_match = min((abs(r - target), r) for r in ratios) if ratios else (float('inf'), 0)
        # Check if any spacing ratio is close
        spacing_ratios = [abs(spacings[i+1] / spacings[i]) if abs(spacings[i]) > 1e-10 else 0
                         for i in range(len(spacings) - 1)]
        best_spacing_match = min((abs(r - target), r) for r in spacing_ratios) if spacing_ratios else (float('inf'), 0)

        matches[target] = {
            'meaning': meaning,
            'best_ratio_distance': best_ratio_match[0],
            'closest_ratio': best_ratio_match[1],
            'best_spacing_ratio_distance': best_spacing_match[0],
            'closest_spacing_ratio': best_spacing_match[1],
        }

    return {
        'eigenvalues': sorted_eigs,
        'real_parts': real_parts,
        'spacings': spacings,
        'consecutive_ratios': ratios,
        'target_matches': matches,
    }


def scan_coupling_bifurcation(
    primes: List[int],
    gamma_range: Tuple[float, float] = (0.01, 5.0),
    n_gamma: int = 100,
    n_steps: int = 300,
    dt: float = 0.05,
    r_stable: float = 5.0,
    lam: float = 0.05,
) -> Dict[str, object]:
    """
    Scan coupling strength γ to detect bifurcation in collapse dynamics.

    For each γ, run the collapse simulation and record the steady-state
    behavior of ⟨R̂⟩. If bifurcation occurs, the late-time values of ⟨R̂⟩
    will split from 1 attractor to 2, then 4, etc. — following the
    period-doubling route to chaos.

    Returns:
        Dict with 'gammas', 'steady_state_values', 'bifurcation_points'
    """
    gammas = np.linspace(gamma_range[0], gamma_range[1], n_gamma)
    all_steady_states = []

    for gamma in gammas:
        result = run_collapse_simulation(
            primes, gamma=gamma, lam=lam, r_stable=r_stable,
            dt=dt, n_steps=n_steps, use_legendre=True
        )
        # Take late-time resonance values as steady state indicators
        late_resonance = result['resonance_history'][-50:]
        # Record unique values (up to tolerance)
        unique_vals = np.unique(np.round(late_resonance, 3))
        all_steady_states.append({
            'gamma': float(gamma),
            'n_attractors': len(unique_vals),
            'attractor_values': unique_vals.tolist(),
            'mean_resonance': float(np.mean(late_resonance)),
            'std_resonance': float(np.std(late_resonance)),
        })

    # Detect bifurcation points (where n_attractors changes)
    bifurcation_points = []
    for i in range(1, len(all_steady_states)):
        if all_steady_states[i]['n_attractors'] != all_steady_states[i - 1]['n_attractors']:
            bifurcation_points.append({
                'gamma': all_steady_states[i]['gamma'],
                'from_attractors': all_steady_states[i - 1]['n_attractors'],
                'to_attractors': all_steady_states[i]['n_attractors'],
            })

    return {
        'gammas': gammas.tolist(),
        'steady_states': all_steady_states,
        'bifurcation_points': bifurcation_points,
    }


def compute_shell_structure(
    primes: List[int],
    gamma: float = 0.5,
    lam: float = 0.1,
    n_shells: int = 5
) -> Dict[str, object]:
    """
    Simulate collapse at multiple target shells and analyze probability distributions.

    The Symbolic Resonance Atomic Model predicts that different r_stable values
    create discrete "shells" where probability concentrates on specific primes.
    This is analogous to electron shells in hydrogen.

    Returns:
        Dict mapping shell index to {target, dominant_prime, final_probs, entropy}
    """
    shells = {}

    # Target shells centered on the first few primes
    target_values = primes[:n_shells]

    for shell_idx, r_target in enumerate(target_values):
        result = run_collapse_simulation(
            primes, gamma=gamma, lam=lam, r_stable=float(r_target),
            dt=0.05, n_steps=300, use_legendre=True
        )
        shells[shell_idx] = {
            'target': r_target,
            'dominant_prime': result['dominant_prime'],
            'dominant_probability': result['dominant_probability'],
            'final_probs': result['final_probs'].tolist(),
            'final_entropy': result['final_entropy'],
            'collapsed': result['dominant_probability'] > 0.5,
        }

    return shells
