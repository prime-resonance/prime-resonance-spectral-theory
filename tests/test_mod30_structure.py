"""
Tests for the mod-30 sieve structure predictions.

Validates:
1. All primes > 5 fall into exactly the 8 coprime residue classes mod 30
2. The 8 residues partition cleanly into mod-4 classes {Q1, Q3}
3. Prime density is approximately uniform across residue classes (Dirichlet)
4. The mod-90 digital root matrix has the predicted structure
"""

import pytest
from src.prime_utils import (
    sieve_primes, mod30_residue, digital_root, is_coprime_to_30,
    COPRIME_RESIDUES_MOD30, Q3_RESIDUES, Q1_RESIDUES,
)
from src.mod30_structure import (
    classify_primes_by_mod30,
    get_mod4_partition,
    verify_mod30_mod4_alignment,
    generate_coprime_sequence,
    verify_period_24_cycle,
    mod90_digital_root_matrix,
    verify_lateral_90_sums,
    prime_density_by_residue_class,
)


class TestPrimesInMod30Classes:
    """Every prime > 5 must belong to exactly one of the 8 coprime residues mod 30."""

    def test_all_primes_up_to_10000_in_coprime_classes(self):
        """PREDICTION: No prime > 5 has a mod-30 residue outside {1,7,11,13,17,19,23,29}."""
        primes = sieve_primes(10_000)
        for p in primes:
            if p <= 5:
                continue
            r = mod30_residue(p)
            assert r in COPRIME_RESIDUES_MOD30, (
                f"Prime {p} has mod-30 residue {r} which is NOT in the coprime set"
            )

    def test_all_8_classes_populated(self):
        """PREDICTION: Each of the 8 coprime residue classes contains primes."""
        primes = sieve_primes(1000)
        classes = classify_primes_by_mod30(primes)
        for r in sorted(COPRIME_RESIDUES_MOD30):
            assert len(classes[r]) > 0, f"Residue class {r} has no primes (unexpected)"

    def test_coprime_set_exactly_8_elements(self):
        """The coprime residue set has exactly phi(30) = 8 elements."""
        assert len(COPRIME_RESIDUES_MOD30) == 8

    def test_coprime_pairs_sum_to_30(self):
        """
        PREDICTION: The 8 coprime residues form 4 complementary pairs summing to 30.
        From the formalism: 1+29=30, 7+23=30, 11+19=30, 13+17=30
        """
        expected_pairs = [(1, 29), (7, 23), (11, 19), (13, 17)]
        for a, b in expected_pairs:
            assert a + b == 30, f"Pair ({a}, {b}) doesn't sum to 30"
            assert a in COPRIME_RESIDUES_MOD30
            assert b in COPRIME_RESIDUES_MOD30


class TestMod4Partition:
    """The mod-30 residues partition into two mod-4 groups governing Legendre asymmetry."""

    def test_q3_residues_all_congruent_3_mod4(self):
        """PREDICTION: {7, 11, 19, 23} ≡ 3 (mod 4)."""
        for r in Q3_RESIDUES:
            assert r % 4 == 3, f"Residue {r} is NOT ≡ 3 (mod 4)"

    def test_q1_residues_all_congruent_1_mod4(self):
        """PREDICTION: {1, 13, 17, 29} ≡ 1 (mod 4)."""
        for r in Q1_RESIDUES:
            assert r % 4 == 1, f"Residue {r} is NOT ≡ 1 (mod 4)"

    def test_partition_is_complete(self):
        """Q1 ∪ Q3 = full coprime residue set."""
        assert Q1_RESIDUES | Q3_RESIDUES == COPRIME_RESIDUES_MOD30

    def test_partition_is_disjoint(self):
        """Q1 ∩ Q3 = ∅."""
        assert len(Q1_RESIDUES & Q3_RESIDUES) == 0

    def test_partition_is_balanced(self):
        """Each partition has exactly 4 elements."""
        assert len(Q1_RESIDUES) == 4
        assert len(Q3_RESIDUES) == 4

    def test_actual_primes_mod60_determines_mod4(self):
        """
        REFINED PREDICTION: Mod-30 class alone does NOT determine mod-4 class.
        E.g., 31 ≡ 1 (mod 30) but 31 ≡ 3 (mod 4); 61 ≡ 1 (mod 30) and 61 ≡ 1 (mod 4).

        The CORRECT partition requires mod-60 (= LCM(30, 4)):
        - Primes ≡ r (mod 60) where r ≡ 3 (mod 4) → Q3
        - Primes ≡ r (mod 60) where r ≡ 1 (mod 4) → Q1

        Within each mod-30 class, approximately half the primes are Q3 and half Q1.
        """
        primes = sieve_primes(5000)
        # Count Q1 vs Q3 primes within each mod-30 class
        q3_by_r30 = {r: 0 for r in COPRIME_RESIDUES_MOD30}
        q1_by_r30 = {r: 0 for r in COPRIME_RESIDUES_MOD30}

        for p in primes:
            if p <= 5:
                continue
            r30 = mod30_residue(p)
            if r30 not in COPRIME_RESIDUES_MOD30:
                continue
            if p % 4 == 3:
                q3_by_r30[r30] += 1
            else:
                q1_by_r30[r30] += 1

        # Each mod-30 class should contain BOTH Q1 and Q3 primes
        for r in sorted(COPRIME_RESIDUES_MOD30):
            total = q3_by_r30[r] + q1_by_r30[r]
            if total > 10:  # Need enough primes for meaningful check
                assert q3_by_r30[r] > 0, f"Mod-30 class {r} has no Q3 primes"
                assert q1_by_r30[r] > 0, f"Mod-30 class {r} has no Q1 primes"
                # Approximately balanced (within 20% of 50/50)
                fraction_q3 = q3_by_r30[r] / total
                assert 0.3 < fraction_q3 < 0.7, (
                    f"Mod-30 class {r}: Q3 fraction = {fraction_q3:.3f}, "
                    f"expected ~0.5"
                )


class TestPrimeDensity:
    """Dirichlet's theorem: primes are approximately equidistributed across coprime classes."""

    def test_density_approximately_uniform(self):
        """
        PREDICTION: Each of the 8 classes contains approximately 12.5% of primes.
        We allow a 3% tolerance for finite-size effects.
        """
        densities = prime_density_by_residue_class(100_000)
        expected = 1.0 / 8  # = 0.125
        tolerance = 0.03

        for r, density in densities.items():
            assert abs(density - expected) < tolerance, (
                f"Residue class {r} has density {density:.4f}, "
                f"expected {expected:.4f} ± {tolerance}"
            )


class TestMod90Matrix:
    """The mod-90 digital root matrix has specific structural properties."""

    def test_matrix_has_6_digital_root_groups(self):
        """PREDICTION: Only digital roots {1, 2, 4, 5, 7, 8} appear."""
        matrix = mod90_digital_root_matrix()
        assert set(matrix.keys()) == {1, 2, 4, 5, 7, 8}

    def test_each_group_has_4_residues(self):
        """PREDICTION: Each digital root maps to exactly 4 mod-90 residues."""
        matrix = mod90_digital_root_matrix()
        for dr, residues in matrix.items():
            assert len(residues) == 4, (
                f"Digital root {dr} has {len(residues)} residues, expected 4"
            )

    def test_lateral_90_sums(self):
        """
        PREDICTION: The 24 mod-90 residues form 12 pairs summing to 90.
        12 × 90 = 1080 = 360 × 3 (three full rotations of the spiral sieve).
        """
        matrix = mod90_digital_root_matrix()
        sums = verify_lateral_90_sums(matrix)
        assert len(sums) == 12, f"Found {len(sums)} pairs summing to 90, expected 12"

    def test_digital_roots_match_coprime_sequence(self):
        """
        PREDICTION: The first 24 numbers coprime to 30 have digital roots
        matching the period-24 cycle from the formalism.
        """
        seq = generate_coprime_sequence(24)
        expected_drs = [1, 7, 2, 4, 8, 1, 5, 2, 4, 1, 5, 7, 2, 4, 8, 5, 7, 4, 8, 1, 5, 7, 2, 8]
        actual_drs = [digital_root(n) for n in seq]
        assert actual_drs == expected_drs, (
            f"Digital roots don't match.\nExpected: {expected_drs}\nActual:   {actual_drs}"
        )


class TestPeriod24Cycle:
    """The coprime-to-30 sequence exhibits exact period-24 digital root cycling."""

    def test_period_24_holds_for_3_periods(self):
        """PREDICTION: Digital root pattern repeats exactly every 24 terms for at least 3 periods."""
        is_periodic, roots = verify_period_24_cycle(n_periods=3)
        assert is_periodic, "Period-24 cycle failed to repeat over 3 periods"

    def test_period_24_holds_for_10_periods(self):
        """PREDICTION: Period-24 cycling holds for at least 10 periods (240 terms)."""
        is_periodic, roots = verify_period_24_cycle(n_periods=10)
        assert is_periodic, "Period-24 cycle failed to repeat over 10 periods"

    def test_digital_root_sum_equals_108(self):
        """
        PREDICTION: The sum of one period-24 cycle's digital roots = 108.
        This connects the period-24 structure to the twist unit 108 = 2²×3³.
        """
        from src.prime_utils import PERIOD_24_DIGITAL_ROOTS
        cycle_sum = sum(PERIOD_24_DIGITAL_ROOTS)
        assert cycle_sum == 108, (
            f"Sum of period-24 digital roots = {cycle_sum}, expected 108 = 2²×3³"
        )
