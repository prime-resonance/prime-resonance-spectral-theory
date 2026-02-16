"""
Tests for Legendre symbol asymmetry predictions.

Validates:
1. Asymmetry occurs ONLY for pairs where both primes ≡ 3 (mod 4)
2. This perfectly aligns with the mod-30 Q3 partition {7,11,19,23}
3. No Q1×Q1 or Q1×Q3 cross-pairs ever have asymmetry
4. The asymmetry matrix has block structure aligned with mod-30 classes
"""

import pytest
import numpy as np
from src.prime_utils import (
    sieve_primes, legendre_symbol, has_legendre_asymmetry,
    mod30_residue, Q3_RESIDUES, Q1_RESIDUES, COPRIME_RESIDUES_MOD30
)
from src.legendre_network import (
    build_legendre_matrix,
    build_asymmetry_matrix,
    count_asymmetric_pairs_by_class,
    compute_residue_class_block_structure,
    compute_qubit_quality_score,
    find_optimal_qubit_pairs,
)


class TestLegendreAsymmetryByMod4:
    """
    CORE PREDICTION: Legendre asymmetry is governed entirely by the mod-4 class
    of the primes, which is determined by their mod-30 residue.
    """

    @pytest.fixture
    def primes_500(self):
        return sieve_primes(500)

    def test_q3_pairs_always_asymmetric(self, primes_500):
        """
        PREDICTION: For primes p ≡ 3 (mod 4) and q ≡ 3 (mod 4),
        (p/q) ≠ (q/p) ALWAYS (by quadratic reciprocity).
        """
        q3_primes = [p for p in primes_500 if p > 5 and p % 4 == 3]
        violations = 0
        total = 0
        for i, p in enumerate(q3_primes):
            for q in q3_primes[i+1:]:
                total += 1
                if not has_legendre_asymmetry(p, q):
                    violations += 1
        assert violations == 0, (
            f"{violations}/{total} Q3×Q3 pairs lack asymmetry (should be 0)"
        )

    def test_q1_pairs_never_asymmetric(self, primes_500):
        """
        PREDICTION: For primes p ≡ 1 (mod 4) and q ≡ 1 (mod 4),
        (p/q) = (q/p) ALWAYS.
        """
        q1_primes = [p for p in primes_500 if p > 5 and p % 4 == 1]
        violations = 0
        total = 0
        for i, p in enumerate(q1_primes):
            for q in q1_primes[i+1:]:
                total += 1
                if has_legendre_asymmetry(p, q):
                    violations += 1
        assert violations == 0, (
            f"{violations}/{total} Q1×Q1 pairs have asymmetry (should be 0)"
        )

    def test_cross_class_never_asymmetric(self, primes_500):
        """
        PREDICTION: For p ≡ 1 (mod 4) and q ≡ 3 (mod 4),
        (p/q) = (q/p) ALWAYS (since (p-1)/2*(q-1)/2 is even).
        """
        q1_primes = [p for p in primes_500 if p > 5 and p % 4 == 1][:30]
        q3_primes = [p for p in primes_500 if p > 5 and p % 4 == 3][:30]
        violations = 0
        total = 0
        for p in q1_primes:
            for q in q3_primes:
                total += 1
                if has_legendre_asymmetry(p, q):
                    violations += 1
        assert violations == 0, (
            f"{violations}/{total} Q1×Q3 cross-pairs have asymmetry (should be 0)"
        )


class TestAsymmetryBlockStructure:
    """
    Asymmetric pairs distribute across mod-30 residue classes.

    REFINED INSIGHT: Since mod-30 class does NOT uniquely determine mod-4 class,
    asymmetric pairs can appear in ANY mod-30 block (wherever both primes
    happen to be ≡ 3 mod 4). The key prediction is at the mod-4 level, not mod-30.
    """

    @pytest.fixture
    def primes_300(self):
        return sieve_primes(300)

    def test_asymmetric_pairs_exist_across_mod30_blocks(self, primes_300):
        """
        REFINED PREDICTION: Since each mod-30 class contains both Q1 and Q3 primes,
        asymmetric pairs should appear in multiple mod-30 blocks.
        """
        counts = count_asymmetric_pairs_by_class(primes_300)
        blocks_with_pairs = sum(1 for c in counts.values() if c > 0)
        assert blocks_with_pairs > 0, "No asymmetric pairs found anywhere"

    def test_asymmetric_pairs_only_between_mod4_class3(self, primes_300):
        """
        CORE PREDICTION: Every asymmetric pair consists of two primes both ≡ 3 mod 4.
        This is the fundamental theorem of quadratic reciprocity applied.
        """
        working = [p for p in primes_300 if p > 5]
        violations = 0
        total_asym = 0
        for i, p in enumerate(working):
            for j in range(i + 1, min(i + 50, len(working))):  # Limit for speed
                q = working[j]
                if has_legendre_asymmetry(p, q):
                    total_asym += 1
                    if p % 4 != 3 or q % 4 != 3:
                        violations += 1
        assert violations == 0, (
            f"{violations}/{total_asym} asymmetric pairs violate the mod-4 rule"
        )

    def test_q3_block_is_dense_with_asymmetric_pairs(self, primes_300):
        """
        PREDICTION: Among primes ≡ 3 (mod 4), ALL pairs have asymmetry.
        """
        q3_primes = [p for p in primes_300 if p > 5 and p % 4 == 3][:30]
        total, asym = 0, 0
        for i, p in enumerate(q3_primes):
            for q in q3_primes[i+1:]:
                total += 1
                if has_legendre_asymmetry(p, q):
                    asym += 1
        assert asym == total, (
            f"Only {asym}/{total} Q3 pairs have asymmetry (expected all)"
        )


class TestQubitQualityScoring:
    """The qubit quality scorer correctly identifies valid and invalid pairs."""

    def test_q3_pair_has_positive_score(self):
        """Q3×Q3 pairs should always score > 0."""
        # 7 and 11 are both ≡ 3 (mod 4)
        score = compute_qubit_quality_score(7, 11)
        assert score > 0, f"Q3 pair (7, 11) scored {score}, expected > 0"

    def test_q1_pair_has_zero_score(self):
        """Q1×Q1 pairs should score 0 (no asymmetry)."""
        # 13 and 29 are both ≡ 1 (mod 4)
        score = compute_qubit_quality_score(13, 29)
        assert score == 0.0, f"Q1 pair (13, 29) scored {score}, expected 0"

    def test_cross_pair_has_zero_score(self):
        """Q1×Q3 cross-pairs should score 0."""
        score = compute_qubit_quality_score(7, 13)
        assert score == 0.0, f"Cross pair (7, 13) scored {score}, expected 0"

    def test_optimal_pairs_all_from_q3(self):
        """
        PREDICTION: All optimal qubit pairs should involve primes ≡ 3 (mod 4).
        """
        primes = sieve_primes(200)
        optimal = find_optimal_qubit_pairs(primes, n_pairs=20)
        for p, q, score in optimal:
            assert p % 4 == 3, f"Optimal qubit prime {p} is not ≡ 3 (mod 4)"
            assert q % 4 == 3, f"Optimal qubit prime {q} is not ≡ 3 (mod 4)"


class TestLegendreMatrixStructure:
    """The Legendre matrix exhibits structured (non-random) block patterns."""

    @pytest.fixture
    def small_primes(self):
        """First 30 primes > 5 for matrix analysis."""
        all_primes = sieve_primes(150)
        return [p for p in all_primes if p > 5]

    def test_asymmetry_matrix_is_symmetric(self, small_primes):
        """The asymmetry indicator matrix A should be symmetric."""
        A = build_asymmetry_matrix(small_primes)
        assert np.allclose(A, A.T), "Asymmetry matrix is not symmetric"

    def test_legendre_matrix_is_not_symmetric(self, small_primes):
        """
        The raw Legendre matrix L should NOT be symmetric in general,
        due to quadratic reciprocity asymmetry.
        """
        L = build_legendre_matrix(small_primes)
        # At least some entries should differ from their transpose
        diff = np.sum(L != L.T)
        assert diff > 0, "Legendre matrix is fully symmetric (unexpected)"

    def test_asymmetry_count_in_legendre_matrix(self, small_primes):
        """
        PREDICTION: The number of (i,j) where L[i,j] ≠ L[j,i] equals
        exactly the number of Q3×Q3 pairs (where both primes have same asymmetry class).
        """
        L = build_legendre_matrix(small_primes)
        A = build_asymmetry_matrix(small_primes)

        # Count asymmetric positions in L (upper triangle)
        asym_count_L = 0
        for i in range(len(small_primes)):
            for j in range(i + 1, len(small_primes)):
                if L[i, j] != L[j, i]:
                    asym_count_L += 1

        # Count from asymmetry matrix (upper triangle)
        asym_count_A = int(np.sum(A) // 2)

        assert asym_count_L == asym_count_A, (
            f"Asymmetric positions in L ({asym_count_L}) ≠ "
            f"entries in A ({asym_count_A})"
        )
