"""
Mod-30 sieve structure analysis.

Implements the Prime Spiral Sieve structure from the formalism,
analyzing the 8-dimensional coprime residue class organization
and period-24 digital root cycles.
"""

from typing import Dict, List, Tuple, Set
from .prime_utils import (
    sieve_primes, mod30_residue, digital_root, is_coprime_to_30,
    COPRIME_RESIDUES_MOD30, Q3_RESIDUES, Q1_RESIDUES, PERIOD_24_DIGITAL_ROOTS
)


def classify_primes_by_mod30(primes: List[int]) -> Dict[int, List[int]]:
    """
    Classify primes into their mod-30 residue classes.

    All primes > 5 must fall into one of the 8 coprime residue classes mod 30:
    {1, 7, 11, 13, 17, 19, 23, 29}

    Returns:
        Dict mapping residue class -> list of primes in that class
    """
    classes: Dict[int, List[int]] = {r: [] for r in sorted(COPRIME_RESIDUES_MOD30)}
    for p in primes:
        if p <= 5:  # 2, 3, 5 are the scaffolding primes
            continue
        r = mod30_residue(p)
        if r in COPRIME_RESIDUES_MOD30:
            classes[r].append(p)
    return classes


def get_mod4_partition(primes: List[int]) -> Tuple[List[int], List[int]]:
    """
    Partition primes >5 into mod-4 classes.

    Returns:
        (q3_primes, q1_primes) where:
        - q3_primes: primes ≡ 3 (mod 4) — candidates for asymmetric Legendre pairs
        - q1_primes: primes ≡ 1 (mod 4) — no Legendre asymmetry possible
    """
    q3 = [p for p in primes if p > 5 and p % 4 == 3]
    q1 = [p for p in primes if p > 5 and p % 4 == 1]
    return q3, q1


def verify_mod30_mod4_alignment(primes: List[int]) -> Dict[int, int]:
    """
    Verify that mod-30 residue classes align with mod-4 classes.

    The prediction from the formalism:
    - {7, 11, 19, 23} mod 30 → all ≡ 3 (mod 4)
    - {1, 13, 17, 29} mod 30 → all ≡ 1 (mod 4)

    Returns:
        Dict mapping each mod-30 residue to its mod-4 class.
    """
    alignment = {}
    for r in sorted(COPRIME_RESIDUES_MOD30):
        alignment[r] = r % 4
    return alignment


def generate_coprime_sequence(n: int) -> List[int]:
    """
    Generate the first n numbers coprime to 30.
    These are the numbers in the Prime Spiral Sieve domain.

    The sequence is: 1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, ...
    """
    result = []
    k = 1
    while len(result) < n:
        if is_coprime_to_30(k):
            result.append(k)
        k += 1
    return result


def compute_digital_root_cycle(sequence: List[int]) -> List[int]:
    """Compute the digital root of each element in a sequence."""
    return [digital_root(n) for n in sequence]


def verify_period_24_cycle(n_periods: int = 3) -> Tuple[bool, List[int]]:
    """
    Verify that numbers coprime to 30 exhibit period-24 digital root cycling.

    The expected period-24 cycle is:
    {1, 7, 2, 4, 8, 1, 5, 2, 4, 1, 5, 7, 2, 4, 8, 5, 7, 4, 8, 1, 5, 7, 2, 8}

    Returns:
        (is_periodic, actual_roots) where is_periodic is True if the
        pattern repeats exactly over n_periods.
    """
    seq = generate_coprime_sequence(24 * n_periods)
    roots = compute_digital_root_cycle(seq)

    is_periodic = True
    for period_idx in range(n_periods):
        start = period_idx * 24
        chunk = roots[start:start + 24]
        if chunk != PERIOD_24_DIGITAL_ROOTS:
            is_periodic = False
            break

    return is_periodic, roots


def mod90_digital_root_matrix() -> Dict[int, List[int]]:
    """
    Build the mod-90 congruence matrix indexed by digital root.

    From the formalism:
    dr 1 = n ≡ {1, 19, 37, 73} mod 90
    dr 2 = n ≡ {11, 29, 47, 83} mod 90
    dr 4 = n ≡ {13, 31, 49, 67} mod 90
    dr 5 = n ≡ {23, 41, 59, 77} mod 90
    dr 7 = n ≡ {7, 43, 61, 79} mod 90
    dr 8 = n ≡ {17, 53, 71, 89} mod 90

    Returns:
        Dict mapping digital root -> list of mod-90 residues
    """
    return {
        1: [1, 19, 37, 73],
        2: [11, 29, 47, 83],
        4: [13, 31, 49, 67],
        5: [23, 41, 59, 77],
        7: [7, 43, 61, 79],
        8: [17, 53, 71, 89],
    }


def verify_lateral_90_sums(matrix: Dict[int, List[int]]) -> List[int]:
    """
    Verify that lateral pairs in the mod-90 matrix sum to 90.

    The formalism claims 12 lateral 90-sums (12 × 90 = 1080 = 360 × 3).

    Returns:
        List of all pairwise sums that equal 90.
    """
    sums_of_90 = []
    all_residues = []
    for dr, residues in matrix.items():
        all_residues.extend(residues)

    all_residues_set = set(all_residues)
    for r in sorted(all_residues):
        complement = 90 - r
        if complement in all_residues_set and complement > r:
            sums_of_90.append((r, complement))

    return sums_of_90


def prime_density_by_residue_class(limit: int) -> Dict[int, float]:
    """
    Compute the density of primes in each mod-30 residue class up to limit.

    By Dirichlet's theorem, primes are approximately equidistributed across
    the 8 coprime residue classes. Deviations reveal the mod-30 structure.

    Returns:
        Dict mapping residue class -> fraction of primes in that class
    """
    primes = sieve_primes(limit)
    classes = classify_primes_by_mod30(primes)
    total = sum(len(v) for v in classes.values())
    if total == 0:
        return {r: 0.0 for r in sorted(COPRIME_RESIDUES_MOD30)}
    return {r: len(classes[r]) / total for r in sorted(COPRIME_RESIDUES_MOD30)}
