"""
Spectral Constants: Emergence of physical ratios from the Legendre Hamiltonian.

The twist framework claims all physical constants emerge from f(108, trefoil, primes).
This module tests whether the eigenvalue spectrum of the Legendre-weighted Hamiltonian,
as the prime basis grows, exhibits scaling laws and ratios that converge toward
physical constants.

Key predictions tested:
1. Spectral scaling: how eigenvalue spread grows with basis size
2. Ground state energy relates to log(p) structure
3. Eigenvalue spacing statistics deviate from random matrices (GUE)
4. The ratio λ_max/λ_min approaches values related to 108 or 17
5. The spectral density shows structure at positions related to 30, 108, 288
"""

from typing import List, Tuple, Dict
import numpy as np
from scipy.linalg import eigvalsh
from .prime_utils import sieve_primes, digital_root, TWIST_UNIT
from .collapse_hamiltonian import (
    build_full_hamiltonian,
    build_kinetic_operator,
    build_resonance_potential,
    build_legendre_resonance_potential,
)


def spectral_scaling_law(
    max_prime: int = 500,
    n_samples: int = 10,
    gamma: float = 1.0,
    use_legendre: bool = True,
) -> Dict[str, object]:
    """
    Study how the eigenvalue spectrum scales as we increase the prime basis size.

    For each sample size N, we take the first N primes and compute the
    Hamiltonian's eigenvalue spectrum. We track:
    - Spectral width (max - min real eigenvalue)
    - Spectral centroid (mean eigenvalue)
    - Standard deviation of eigenvalues
    - Ground state energy (most negative eigenvalue)

    Returns:
        Dict with 'basis_sizes', 'widths', 'centroids', 'stds', 'ground_states'
    """
    all_primes = sieve_primes(max_prime)
    all_primes = [p for p in all_primes if p > 2]  # Exclude p=2 for cleaner Legendre

    # Sample at logarithmically spaced sizes
    max_n = len(all_primes)
    sizes = sorted(set(
        [min(max_n, max(4, int(np.exp(k)))) for k in np.linspace(np.log(4), np.log(max_n), n_samples)]
    ))

    results = {
        'basis_sizes': [], 'widths': [], 'centroids': [],
        'stds': [], 'ground_states': [], 'max_eigs': [],
    }

    for n in sizes:
        primes = all_primes[:n]
        H = build_full_hamiltonian(primes, gamma, use_legendre)
        eigs = np.linalg.eigvals(H)
        real_eigs = np.sort(eigs.real)

        results['basis_sizes'].append(n)
        results['widths'].append(float(real_eigs[-1] - real_eigs[0]))
        results['centroids'].append(float(np.mean(real_eigs)))
        results['stds'].append(float(np.std(real_eigs)))
        results['ground_states'].append(float(real_eigs[0]))
        results['max_eigs'].append(float(real_eigs[-1]))

    return results


def eigenvalue_ratio_analysis(
    primes: List[int],
    gamma: float = 1.0,
    use_legendre: bool = True,
) -> Dict[str, object]:
    """
    Analyze eigenvalue ratios looking for physical constant signatures.

    Key ratios to search for:
    - 17 (trefoil complexity, m_p/m_e factor)
    - 29 (boundary prime correction to α)
    - 108 (twist unit)
    - 137 (fine structure constant)

    Returns:
        Dict with eigenvalue data and closest matches to target constants.
    """
    H = build_full_hamiltonian(primes, gamma, use_legendre)
    eigs = np.linalg.eigvals(H)
    real_eigs = np.sort(eigs.real)

    # Compute all pairwise ratios (magnitude)
    n = len(real_eigs)
    all_ratios = []
    for i in range(n):
        for j in range(i + 1, n):
            if abs(real_eigs[i]) > 0.01:
                all_ratios.append(abs(real_eigs[j] / real_eigs[i]))

    # Spacings
    spacings = np.diff(real_eigs)
    spacing_ratios = []
    for i in range(len(spacings) - 1):
        if abs(spacings[i]) > 1e-6:
            spacing_ratios.append(abs(spacings[i + 1] / spacings[i]))

    # Target physical constants
    targets = {
        'trefoil_17': 17,
        'boundary_29': 29,
        'twist_108': 108,
        'alpha_137': 137,
        'width_to_108_ratio': None,  # Will compute
    }

    # Find closest matches
    matches = {}
    for name, target in targets.items():
        if target is None:
            continue
        if all_ratios:
            diffs = [abs(r - target) for r in all_ratios]
            best_idx = np.argmin(diffs)
            matches[name] = {
                'target': target,
                'closest_ratio': all_ratios[best_idx],
                'distance': diffs[best_idx],
                'relative_error': diffs[best_idx] / target,
            }

    # Width-to-108 analysis
    width = real_eigs[-1] - real_eigs[0]
    matches['width_mod_108'] = {
        'width': float(width),
        'width_div_108': float(width / TWIST_UNIT) if TWIST_UNIT != 0 else 0,
        'width_mod_108': float(width % TWIST_UNIT),
    }

    return {
        'eigenvalues': real_eigs.tolist(),
        'n_eigenvalues': n,
        'spectral_width': float(real_eigs[-1] - real_eigs[0]),
        'all_ratios': all_ratios,
        'spacing_ratios': spacing_ratios,
        'matches': matches,
    }


def wigner_surmise_comparison(
    primes: List[int],
    gamma: float = 1.0,
    use_legendre: bool = True,
) -> Dict[str, float]:
    """
    Compare eigenvalue spacing distribution to the Wigner surmise (GUE prediction).

    For random matrices (GUE), the nearest-neighbor spacing distribution follows:
        P(s) = (32/π²) s² exp(-4s²/π)

    Deviation from this indicates non-random structure in the Hamiltonian.

    Returns:
        Dict with 'mean_spacing', 'std_spacing', 'skewness',
        'gue_expected_mean', 'gue_expected_variance',
        'deviation_from_gue' (chi-squared statistic)
    """
    H = build_full_hamiltonian(primes, gamma, use_legendre)
    eigs = np.linalg.eigvals(H)
    real_eigs = np.sort(eigs.real)

    spacings = np.diff(real_eigs)
    if len(spacings) == 0:
        return {'mean_spacing': 0, 'std_spacing': 0, 'deviation_from_gue': 0}

    mean_s = np.mean(spacings)
    if mean_s == 0:
        return {'mean_spacing': 0, 'std_spacing': 0, 'deviation_from_gue': 0}

    # Normalize spacings
    norm_spacings = spacings / mean_s

    # GUE statistics
    gue_mean = 1.0  # By normalization
    gue_var = (4 - np.pi) / np.pi  # ≈ 0.273 for GUE

    observed_mean = float(np.mean(norm_spacings))
    observed_var = float(np.var(norm_spacings))

    # Simple deviation measure: how far is the observed variance from GUE?
    variance_deviation = abs(observed_var - gue_var) / gue_var

    # Level repulsion: GUE has P(0) = 0 (level repulsion).
    # Count near-zero spacings
    n_small = np.sum(norm_spacings < 0.1)
    small_fraction = n_small / len(norm_spacings)

    # GUE predicts very few s ≈ 0 spacings
    gue_small_fraction = 0.01  # Approximate
    repulsion_deviation = abs(small_fraction - gue_small_fraction)

    return {
        'mean_spacing': float(mean_s),
        'std_spacing': float(np.std(norm_spacings)),
        'observed_variance': observed_var,
        'gue_expected_variance': gue_var,
        'variance_deviation': variance_deviation,
        'small_spacing_fraction': float(small_fraction),
        'repulsion_deviation': repulsion_deviation,
    }


def trace_invariants(
    primes: List[int],
    gamma: float = 1.0,
    use_legendre: bool = True,
    max_power: int = 6,
) -> Dict[int, complex]:
    """
    Compute trace invariants Tr(H^k) for k = 1, ..., max_power.

    These are spectral invariants related to closed walks in the interaction graph.
    The formalism's 288-element symmetry predicts specific relationships:
    - Tr(H^2) relates to total coupling strength
    - Tr(H^3) detects triangular (3-body) correlations
    - The ratio Tr(H^3)/Tr(H^2) may relate to the trefoil crossing number c=3

    Returns:
        Dict mapping power k -> Tr(H^k)
    """
    H = build_full_hamiltonian(primes, gamma, use_legendre)
    traces = {}
    H_power = np.eye(len(primes), dtype=np.complex128)

    for k in range(1, max_power + 1):
        H_power = H_power @ H
        traces[k] = complex(np.trace(H_power))

    return traces


def mod30_block_spectrum(
    primes: List[int],
    gamma: float = 1.0,
) -> Dict[int, List[float]]:
    """
    Compute the spectrum of Hamiltonian blocks restricted to each mod-30 class.

    If the mod-30 structure organizes the Hamiltonian into near-independent blocks,
    each block's eigenvalues should differ systematically.

    Returns:
        Dict mapping mod-30 residue -> eigenvalues of the restricted block
    """
    from .prime_utils import mod30_residue, COPRIME_RESIDUES_MOD30

    # Group primes by mod-30 class
    groups = {r: [] for r in sorted(COPRIME_RESIDUES_MOD30)}
    prime_to_idx = {p: i for i, p in enumerate(primes)}

    for p in primes:
        r = mod30_residue(p)
        if r in COPRIME_RESIDUES_MOD30:
            groups[r].append(p)

    block_spectra = {}
    for r, group_primes in groups.items():
        if len(group_primes) < 2:
            block_spectra[r] = []
            continue

        # Build sub-Hamiltonian for this block
        n = len(group_primes)
        H_block = np.zeros((n, n), dtype=np.complex128)
        for i in range(n):
            H_block[i, i] = -1j * np.log(group_primes[i])
            for j in range(n):
                if i != j:
                    ls = 1  # Default for same-block coupling
                    try:
                        from .prime_utils import legendre_symbol
                        ls = legendre_symbol(group_primes[i], group_primes[j])
                    except (ValueError, ZeroDivisionError):
                        ls = 1
                    H_block[i, j] = -gamma * np.log(group_primes[i] * group_primes[j]) * ls

        eigs = np.linalg.eigvals(H_block)
        block_spectra[r] = sorted(eigs.real.tolist())

    return block_spectra
