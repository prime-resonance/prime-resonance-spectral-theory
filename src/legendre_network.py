"""
Legendre symbol interaction network.

Constructs the Legendre interaction matrix L(p_i, p_j) = (p_i / p_j) for
a collection of primes and analyzes its graph-theoretic properties.
"""

from typing import Dict, List, Tuple, Optional
import numpy as np
from .prime_utils import (
    sieve_primes, legendre_symbol, has_legendre_asymmetry,
    mod30_residue, COPRIME_RESIDUES_MOD30, Q3_RESIDUES, Q1_RESIDUES
)


def build_legendre_matrix(primes: List[int]) -> np.ndarray:
    """
    Build the Legendre interaction matrix for a list of primes.

    L[i,j] = legendre_symbol(primes[i], primes[j]) for i ≠ j
    L[i,i] = 0 (no self-interaction)

    This matrix is generally NOT symmetric due to quadratic reciprocity:
    L[i,j] * L[j,i] = (-1)^((p_i-1)/2 * (p_j-1)/2)

    Returns:
        np.ndarray of shape (n, n) with entries in {-1, 0, +1}
    """
    n = len(primes)
    L = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i != j:
                L[i, j] = legendre_symbol(primes[i], primes[j])
    return L


def build_asymmetry_matrix(primes: List[int]) -> np.ndarray:
    """
    Build the asymmetry indicator matrix.

    A[i,j] = 1 if (p_i/p_j) ≠ (p_j/p_i), else 0

    This matrix is symmetric. Non-zero entries indicate pairs that can
    form topological qubits (from formalism_v2_logic.md).

    Returns:
        np.ndarray of shape (n, n) with entries in {0, 1}
    """
    n = len(primes)
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if has_legendre_asymmetry(primes[i], primes[j]):
                A[i, j] = 1
                A[j, i] = 1
    return A


def build_resonance_hamiltonian(primes: List[int]) -> np.ndarray:
    """
    Build the resonance Hamiltonian matrix.

    H[i,j] = (p_i/p_j) + (p_j/p_i) for i ≠ j

    When (p_i/p_j) = (p_j/p_i), H[i,j] = ±2 (strong coupling)
    When (p_i/p_j) ≠ (p_j/p_i), H[i,j] = 0 (frustrated/asymmetric)

    This is a symmetric matrix suitable for eigenvalue analysis.

    Returns:
        np.ndarray of shape (n, n)
    """
    n = len(primes)
    H = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            ls_ij = legendre_symbol(primes[i], primes[j])
            ls_ji = legendre_symbol(primes[j], primes[i])
            H[i, j] = ls_ij + ls_ji
            H[j, i] = H[i, j]
    return H


def compute_residue_class_block_structure(
    primes: List[int],
    legendre_mat: np.ndarray
) -> Dict[Tuple[int, int], float]:
    """
    Compute the mean Legendre symbol value for each pair of mod-30 residue classes.

    If primes in the same residue class have correlated Legendre relationships
    with third-party primes, we expect the block means to differ significantly
    from zero and from each other.

    Returns:
        Dict mapping (residue_a, residue_b) -> mean Legendre symbol value
    """
    # Map each prime to its index and residue class
    prime_to_idx = {p: i for i, p in enumerate(primes)}
    residue_groups: Dict[int, List[int]] = {r: [] for r in sorted(COPRIME_RESIDUES_MOD30)}

    for p in primes:
        if p <= 5:
            continue
        r = mod30_residue(p)
        if r in COPRIME_RESIDUES_MOD30:
            residue_groups[r].append(prime_to_idx[p])

    block_means = {}
    residues = sorted(COPRIME_RESIDUES_MOD30)
    for ra in residues:
        for rb in residues:
            indices_a = residue_groups[ra]
            indices_b = residue_groups[rb]
            if not indices_a or not indices_b:
                block_means[(ra, rb)] = 0.0
                continue

            total = 0.0
            count = 0
            for i in indices_a:
                for j in indices_b:
                    if i != j:
                        total += legendre_mat[i, j]
                        count += 1
            block_means[(ra, rb)] = total / count if count > 0 else 0.0

    return block_means


def count_asymmetric_pairs_by_class(
    primes: List[int]
) -> Dict[Tuple[int, int], int]:
    """
    Count asymmetric Legendre pairs grouped by (mod-30 class, mod-30 class).

    The prediction: asymmetric pairs should ONLY appear when both primes
    are in Q3_RESIDUES = {7, 11, 19, 23} (both ≡ 3 mod 4).

    Returns:
        Dict mapping (residue_a, residue_b) -> count of asymmetric pairs
    """
    # Filter to primes > 5
    working_primes = [p for p in primes if p > 5]
    residues = sorted(COPRIME_RESIDUES_MOD30)

    # Group primes by mod-30 residue
    groups: Dict[int, List[int]] = {r: [] for r in residues}
    for p in working_primes:
        r = mod30_residue(p)
        if r in COPRIME_RESIDUES_MOD30:
            groups[r].append(p)

    counts: Dict[Tuple[int, int], int] = {}
    for ra in residues:
        for rb in residues:
            if ra > rb:
                continue
            count = 0
            for p in groups[ra]:
                target_group = groups[rb] if ra != rb else [q for q in groups[rb] if q > p]
                for q in target_group:
                    if has_legendre_asymmetry(p, q):
                        count += 1
            counts[(ra, rb)] = count

    return counts


def compute_qubit_quality_score(p: int, q: int) -> float:
    """
    Compute a quality score for a prime pair as a topological qubit.

    A good qubit pair has:
    1. Legendre asymmetry: (p/q) ≠ (q/p)
    2. Both symbols non-zero
    3. Clear separation of eigenvalues (|p - q| not too small or too large)

    Returns:
        Float score ≥ 0 (higher is better)
    """
    ls_pq = legendre_symbol(p, q)
    ls_qp = legendre_symbol(q, p)

    if ls_pq == 0 or ls_qp == 0:
        return 0.0

    # Asymmetry is essential
    if ls_pq == ls_qp:
        return 0.0

    # Quality factor: prefer moderate separation
    ratio = min(p, q) / max(p, q)
    quality = ratio * (1 - ratio)  # Peaks at ratio = 0.5

    return quality


def find_optimal_qubit_pairs(
    primes: List[int],
    n_pairs: int = 10
) -> List[Tuple[int, int, float]]:
    """
    Find the top-scoring prime pairs for use as topological qubits.

    Returns:
        List of (p, q, score) tuples, sorted by score descending.
    """
    pairs = []
    for i, p in enumerate(primes):
        if p <= 5:
            continue
        for j in range(i + 1, len(primes)):
            q = primes[j]
            if q <= 5:
                continue
            score = compute_qubit_quality_score(p, q)
            if score > 0:
                pairs.append((p, q, score))

    pairs.sort(key=lambda x: x[2], reverse=True)
    return pairs[:n_pairs]
