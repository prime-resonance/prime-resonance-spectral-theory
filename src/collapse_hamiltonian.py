"""
Prime Resonance Collapse Hamiltonian.

Implements the Symbolic Resonance Atomic Model Hamiltonian from the formalism:
    Ĥ = T̂ + V̂_res

where:
    T̂ = -iℏ Σ_p log(p) |p⟩⟨p|  (kinetic/entropy gradient term)
    V̂_res = Σ_{p≠q} -γ log(pq) |p⟩⟨q|  (resonance potential)

NOVEL EXTENSION: We augment V̂_res with Legendre symbol weighting:
    V̂_leg[i,j] = -γ · log(p_i · p_j) · legendre(p_i, p_j)

This creates a Hamiltonian that encodes both the coupling magnitude (log product)
and the chirality (Legendre symbol sign), unifying the atomic model with the
topological prime computing Legendre structure.

Non-Hermitian collapse dynamics:
    d|Ψ⟩/dt = iĤ|Ψ⟩ - λ(R̂ - r_stable)|Ψ⟩

where R̂ = Σ_p p|p⟩⟨p| is the resonance operator and r_stable is the target shell.
"""

from typing import List, Tuple, Dict, Optional
import numpy as np
from scipy.linalg import expm
from .prime_utils import sieve_primes, legendre_symbol, digital_root


def build_kinetic_operator(primes: List[int]) -> np.ndarray:
    """
    Build the kinetic term T̂ = -i Σ_p log(p) |p⟩⟨p|.

    This is a diagonal matrix with -i*log(p) on the diagonal.
    Represents entropic momentum — primes with larger log carry more entropy.

    Returns:
        Complex diagonal matrix of shape (N, N)
    """
    N = len(primes)
    T = np.zeros((N, N), dtype=np.complex128)
    for i, p in enumerate(primes):
        T[i, i] = -1j * np.log(p)
    return T


def build_resonance_potential(primes: List[int], gamma: float = 1.0) -> np.ndarray:
    """
    Build the standard resonance potential V̂_res = -γ Σ_{p≠q} log(pq)|p⟩⟨q|.

    This is the original form from the Symbolic Resonance Atomic Model.

    Returns:
        Real symmetric matrix of shape (N, N)
    """
    N = len(primes)
    V = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            if i != j:
                V[i, j] = -gamma * np.log(primes[i] * primes[j])
    return V


def build_legendre_resonance_potential(
    primes: List[int],
    gamma: float = 1.0
) -> np.ndarray:
    """
    Build the Legendre-weighted resonance potential (NOVEL).

    V̂_leg[i,j] = -γ · log(p_i · p_j) · legendre(p_i, p_j)

    This combines the magnitude (log product) with the chirality sign (Legendre symbol),
    creating asymmetric coupling that encodes quadratic reciprocity structure.

    Returns:
        Complex matrix of shape (N, N) — NOT symmetric due to Legendre asymmetry
    """
    N = len(primes)
    V = np.zeros((N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            if i != j:
                ls = legendre_symbol(primes[i], primes[j])
                V[i, j] = -gamma * np.log(primes[i] * primes[j]) * ls
    return V


def build_full_hamiltonian(
    primes: List[int],
    gamma: float = 1.0,
    use_legendre: bool = True
) -> np.ndarray:
    """
    Build the full Hamiltonian Ĥ = T̂ + V̂.

    Args:
        primes: List of prime basis states
        gamma: Coupling strength for resonance potential
        use_legendre: If True, use Legendre-weighted potential (novel);
                      if False, use standard log-product potential

    Returns:
        Complex matrix of shape (N, N)
    """
    T = build_kinetic_operator(primes)
    if use_legendre:
        V = build_legendre_resonance_potential(primes, gamma)
    else:
        V = build_resonance_potential(primes, gamma)
    return T + V


def build_resonance_operator(primes: List[int]) -> np.ndarray:
    """
    Build the resonance operator R̂ = Σ_p p|p⟩⟨p|.

    This is a diagonal matrix with prime values on the diagonal.
    Used in the collapse dynamics to define attractor shells.

    Returns:
        Real diagonal matrix of shape (N, N)
    """
    N = len(primes)
    R = np.zeros((N, N), dtype=np.complex128)
    for i, p in enumerate(primes):
        R[i, i] = float(p)
    return R


def build_collapse_hamiltonian(
    primes: List[int],
    gamma: float = 1.0,
    lam: float = 0.1,
    r_stable: float = 5.0,
    use_legendre: bool = True
) -> np.ndarray:
    """
    Build the effective non-Hermitian Hamiltonian for collapse dynamics.

    Ĥ_eff = Ĥ - iλ(R̂ - r_stable)

    The non-Hermitian part drives the system toward states where
    ⟨R̂⟩ = r_stable (the target attractor shell).

    Args:
        primes: List of prime basis states
        gamma: Coupling strength
        lam: Dissipation rate (collapse speed)
        r_stable: Target resonance value (attractor shell)
        use_legendre: Use Legendre-weighted potential

    Returns:
        Complex non-Hermitian matrix of shape (N, N)
    """
    H = build_full_hamiltonian(primes, gamma, use_legendre)
    R = build_resonance_operator(primes)
    N = len(primes)
    identity = np.eye(N, dtype=np.complex128)

    H_eff = H - 1j * lam * (R - r_stable * identity)
    return H_eff


def evolve_state(
    psi: np.ndarray,
    H_eff: np.ndarray,
    dt: float,
    n_steps: int,
    normalize: bool = True
) -> List[np.ndarray]:
    """
    Evolve a quantum state under the effective Hamiltonian.

    |Ψ(t+dt)⟩ = exp(-i Ĥ_eff dt) |Ψ(t)⟩

    For non-Hermitian H_eff, the norm is NOT preserved, so we
    normalize after each step (consistent with non-Hermitian QM).

    Args:
        psi: Initial state vector (complex, normalized)
        H_eff: Effective Hamiltonian
        dt: Time step
        n_steps: Number of evolution steps
        normalize: Whether to renormalize after each step

    Returns:
        List of state vectors at each time step
    """
    U = expm(-1j * H_eff * dt)
    trajectory = [psi.copy()]

    for _ in range(n_steps):
        psi = U @ psi
        if normalize:
            norm = np.linalg.norm(psi)
            if norm > 1e-15:
                psi = psi / norm
        trajectory.append(psi.copy())

    return trajectory


def compute_probabilities(psi: np.ndarray) -> np.ndarray:
    """Compute measurement probabilities |⟨p|Ψ⟩|² for each prime basis state."""
    return np.abs(psi) ** 2


def compute_entropy(psi: np.ndarray) -> float:
    """
    Compute the symbolic entropy S(|Ψ⟩) = -Σ |c_p|² log|c_p|².

    Low entropy indicates the state has collapsed to one or few prime modes.
    """
    probs = compute_probabilities(psi)
    probs = probs[probs > 1e-15]  # Avoid log(0)
    return float(-np.sum(probs * np.log(probs)))


def compute_resonance_expectation(psi: np.ndarray, primes: List[int]) -> float:
    """Compute ⟨Ψ|R̂|Ψ⟩ = Σ p |c_p|²."""
    probs = compute_probabilities(psi)
    return float(np.sum(probs * np.array(primes, dtype=np.float64)))


def run_collapse_simulation(
    primes: List[int],
    gamma: float = 0.5,
    lam: float = 0.1,
    r_stable: float = 5.0,
    dt: float = 0.05,
    n_steps: int = 200,
    use_legendre: bool = True
) -> Dict[str, object]:
    """
    Run a complete collapse simulation starting from uniform superposition.

    Args:
        primes: Prime basis states (e.g., [2, 3, 5, 7, 11, 13])
        gamma: Coupling strength
        lam: Dissipation rate
        r_stable: Target attractor shell
        dt: Time step
        n_steps: Number of steps
        use_legendre: Use Legendre-weighted potential

    Returns:
        Dict with 'trajectory', 'entropy_history', 'resonance_history',
        'probability_history', 'final_state', 'dominant_prime'
    """
    N = len(primes)

    # Initial state: uniform superposition
    psi = np.ones(N, dtype=np.complex128) / np.sqrt(N)

    # Build Hamiltonian
    H_eff = build_collapse_hamiltonian(primes, gamma, lam, r_stable, use_legendre)

    # Evolve
    trajectory = evolve_state(psi, H_eff, dt, n_steps)

    # Compute observables
    entropy_history = [compute_entropy(s) for s in trajectory]
    resonance_history = [compute_resonance_expectation(s, primes) for s in trajectory]
    probability_history = [compute_probabilities(s) for s in trajectory]

    # Final state analysis
    final_probs = probability_history[-1]
    dominant_idx = int(np.argmax(final_probs))
    dominant_prime = primes[dominant_idx]

    return {
        'trajectory': trajectory,
        'entropy_history': entropy_history,
        'resonance_history': resonance_history,
        'probability_history': probability_history,
        'final_state': trajectory[-1],
        'final_probs': final_probs,
        'dominant_prime': dominant_prime,
        'dominant_probability': float(final_probs[dominant_idx]),
        'final_entropy': entropy_history[-1],
        'initial_entropy': entropy_history[0],
    }
