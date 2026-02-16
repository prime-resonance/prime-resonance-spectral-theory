"""
Spectral analysis of the Legendre interaction matrix.

Analyzes eigenvalues, spectral gaps, and community structure
to validate predictions about non-random organization.
"""

from typing import Dict, List, Tuple, Optional
import numpy as np
from scipy import linalg
from .prime_utils import sieve_primes, mod30_residue, COPRIME_RESIDUES_MOD30
from .legendre_network import build_legendre_matrix, build_resonance_hamiltonian


def compute_eigenspectrum(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute eigenvalues and eigenvectors of a matrix.

    For asymmetric matrices (Legendre matrix), uses general eigenvalue decomposition.
    For symmetric matrices (resonance Hamiltonian), uses symmetric eigenvalue decomposition.

    Returns:
        (eigenvalues, eigenvectors) — eigenvalues sorted by magnitude
    """
    is_symmetric = np.allclose(matrix, matrix.T)

    if is_symmetric:
        eigenvalues, eigenvectors = linalg.eigh(matrix)
        # Sort by descending magnitude
        idx = np.argsort(-np.abs(eigenvalues))
        return eigenvalues[idx], eigenvectors[:, idx]
    else:
        eigenvalues, eigenvectors = linalg.eig(matrix)
        # Sort by descending magnitude of real part
        idx = np.argsort(-np.abs(eigenvalues.real))
        return eigenvalues[idx], eigenvectors[:, idx]


def spectral_gap(eigenvalues: np.ndarray) -> float:
    """
    Compute the spectral gap: difference between the two largest eigenvalue magnitudes.

    A larger spectral gap indicates stronger community structure and predicts
    faster convergence of resonance-based annealing.

    Returns:
        The spectral gap (non-negative float)
    """
    sorted_mags = np.sort(np.abs(eigenvalues))[::-1]
    if len(sorted_mags) < 2:
        return 0.0
    return float(sorted_mags[0] - sorted_mags[1])


def eigenvalue_spacing_distribution(eigenvalues: np.ndarray) -> np.ndarray:
    """
    Compute the normalized nearest-neighbor spacing distribution of eigenvalues.

    For random matrices (GUE), this follows the Wigner surmise: P(s) = (π/2)s exp(-πs²/4)
    For non-random (structured) matrices, the distribution deviates significantly.

    Returns:
        Array of normalized spacings
    """
    # Use real parts for complex eigenvalues
    real_parts = np.sort(np.real(eigenvalues))

    # Unfolding: normalize spacings by local mean spacing
    spacings = np.diff(real_parts)
    if len(spacings) == 0:
        return np.array([])

    mean_spacing = np.mean(spacings)
    if mean_spacing == 0:
        return spacings

    normalized = spacings / mean_spacing
    return normalized


def modularity_by_residue_class(
    primes: List[int],
    matrix: np.ndarray
) -> float:
    """
    Compute the modularity of the Legendre matrix with respect to mod-30 residue classes.

    Modularity Q measures how much the interaction structure deviates from
    random expectation given the community assignment (mod-30 class).

    Q = (1/2m) Σ_ij [A_ij - k_i*k_j/(2m)] δ(c_i, c_j)

    where A_ij is the adjacency (|L_ij|), k_i is degree, m is total edges,
    and c_i is the community (mod-30 class) of node i.

    Returns:
        Modularity Q ∈ [-0.5, 1.0]. Positive values indicate community structure.
    """
    n = len(primes)
    # Use absolute value of Legendre symbols as adjacency
    A = np.abs(matrix).astype(np.float64)
    np.fill_diagonal(A, 0)

    # Total weight
    m = np.sum(A) / 2
    if m == 0:
        return 0.0

    # Node degrees
    k = np.sum(A, axis=1)

    # Community assignments: mod-30 residue class
    communities = [mod30_residue(p) if p > 5 else 0 for p in primes]

    Q = 0.0
    for i in range(n):
        for j in range(n):
            if communities[i] == communities[j] and communities[i] != 0:
                Q += A[i, j] - k[i] * k[j] / (2 * m)

    Q /= (2 * m)
    return float(Q)


def compare_to_random_matrix(
    primes: List[int],
    n_samples: int = 100
) -> Dict[str, float]:
    """
    Compare the Legendre matrix's spectral properties to random signed matrices.

    Generates n_samples random matrices with the same dimension and entry distribution
    {-1, 0, +1} and compares:
    1. Spectral gap
    2. Largest eigenvalue magnitude
    3. Modularity

    Returns:
        Dict with 'gap_z_score', 'max_eig_z_score', 'modularity_z_score'
        indicating how many standard deviations the real matrix differs from random.
    """
    # Real matrix properties
    L = build_legendre_matrix(primes)
    real_eigenvalues = compute_eigenspectrum(L)[0]
    real_gap = spectral_gap(real_eigenvalues)
    real_max_eig = float(np.max(np.abs(real_eigenvalues)))
    real_modularity = modularity_by_residue_class(primes, L)

    # Random comparison
    n = len(primes)
    rng = np.random.default_rng(42)
    random_gaps = []
    random_max_eigs = []
    random_modularities = []

    for _ in range(n_samples):
        # Random signed matrix with same sparsity pattern
        R = rng.choice([-1, 0, 1], size=(n, n), p=[0.4, 0.2, 0.4])
        np.fill_diagonal(R, 0)

        r_eigs = compute_eigenspectrum(R.astype(np.float64))[0]
        random_gaps.append(spectral_gap(r_eigs))
        random_max_eigs.append(float(np.max(np.abs(r_eigs))))
        random_modularities.append(modularity_by_residue_class(primes, R))

    def z_score(real_val, random_vals):
        mean = np.mean(random_vals)
        std = np.std(random_vals)
        if std == 0:
            return 0.0
        return float((real_val - mean) / std)

    return {
        'gap_z_score': z_score(real_gap, random_gaps),
        'max_eig_z_score': z_score(real_max_eig, random_max_eigs),
        'modularity_z_score': z_score(real_modularity, random_modularities),
        'real_gap': real_gap,
        'real_max_eig': real_max_eig,
        'real_modularity': real_modularity,
    }


def eigenvector_residue_alignment(
    primes: List[int],
    eigenvectors: np.ndarray,
    top_k: int = 4
) -> Dict[int, Dict[int, float]]:
    """
    Analyze how the top eigenvectors align with mod-30 residue classes.

    If the Legendre matrix has community structure aligned with mod-30 classes,
    the leading eigenvectors should show clustering by residue class.

    Returns:
        Dict mapping eigenvector index -> {residue_class: mean_component}
    """
    alignment = {}
    n = len(primes)

    for k in range(min(top_k, eigenvectors.shape[1])):
        vec = np.real(eigenvectors[:, k])
        residue_means: Dict[int, float] = {}
        for r in sorted(COPRIME_RESIDUES_MOD30):
            indices = [i for i, p in enumerate(primes) if p > 5 and mod30_residue(p) == r]
            if indices:
                residue_means[r] = float(np.mean(np.abs(vec[indices])))
            else:
                residue_means[r] = 0.0
        alignment[k] = residue_means

    return alignment
