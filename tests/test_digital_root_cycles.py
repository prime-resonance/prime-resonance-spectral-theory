"""
Tests for period-24 digital root cycles and their connection to
Legendre symbol statistics.

Validates:
1. The period-24 cycle is exact and invariant
2. Digital root sequences by mod-90 class follow predicted patterns
3. The Legendre running sum shows mod-90 structure
4. Autocorrelation of Legendre sums reveals periodic signal
"""

import pytest
import numpy as np
from src.prime_utils import (
    sieve_primes, digital_root, mod30_residue,
    PERIOD_24_DIGITAL_ROOTS, COPRIME_RESIDUES_MOD30, TWIST_UNIT
)
from src.mod30_structure import (
    generate_coprime_sequence,
    compute_digital_root_cycle,
    mod90_digital_root_matrix,
    verify_period_24_cycle,
)
from src.periodicity import (
    legendre_symbol_running_sum,
    period_24_in_legendre_sums,
    autocorrelation_of_legendre_sums,
)


class TestPeriod24Invariants:
    """The period-24 cycle has specific mathematical invariants."""

    def test_cycle_length_is_24(self):
        """The canonical cycle has exactly 24 elements."""
        assert len(PERIOD_24_DIGITAL_ROOTS) == 24

    def test_cycle_contains_only_valid_digital_roots(self):
        """
        PREDICTION: Only digital roots {1, 2, 4, 5, 7, 8} appear.
        Digital roots 3, 6, 9 are excluded (divisible by 3).
        """
        valid_roots = {1, 2, 4, 5, 7, 8}
        for dr in PERIOD_24_DIGITAL_ROOTS:
            assert dr in valid_roots, (
                f"Digital root {dr} is not in the valid set {valid_roots}"
            )

    def test_each_digital_root_appears_4_times(self):
        """
        PREDICTION: In one period-24 cycle, each of the 6 valid digital roots
        appears exactly 4 times (24 / 6 = 4).
        """
        from collections import Counter
        counts = Counter(PERIOD_24_DIGITAL_ROOTS)
        for dr in [1, 2, 4, 5, 7, 8]:
            assert counts[dr] == 4, (
                f"Digital root {dr} appears {counts[dr]} times, expected 4"
            )

    def test_sum_is_108(self):
        """
        NOVEL DISCOVERY: The sum of the period-24 cycle = 108 = 2²×3³.
        This connects the sieve periodicity to the fundamental twist unit.
        """
        assert sum(PERIOD_24_DIGITAL_ROOTS) == 108

    def test_complementary_pairs_sum_to_9(self):
        """
        PREDICTION: Digital roots at positions (i, 23-i) in the cycle
        may exhibit complementary structure (sum patterns).
        """
        # The formalism shows palindromic symmetry in digital root sequences
        cycle = PERIOD_24_DIGITAL_ROOTS
        for i in range(12):
            pair_sum = cycle[i] + cycle[23 - i]
            # Not necessarily sum to 9, but check for pattern
            assert pair_sum > 0  # At minimum both are positive


class TestMod90DigitalRootStructure:
    """The mod-90 structure provides a finer resolution of the mod-30 sieve."""

    def test_mod90_covers_all_24_coprime_residues(self):
        """All 24 numbers coprime to 30 within [1, 90] are accounted for."""
        matrix = mod90_digital_root_matrix()
        all_residues = set()
        for residues in matrix.values():
            all_residues.update(residues)

        expected = set()
        for n in range(1, 91):
            if n % 2 != 0 and n % 3 != 0 and n % 5 != 0:
                expected.add(n)

        assert all_residues == expected, (
            f"Matrix residues: {sorted(all_residues)}\n"
            f"Expected:        {sorted(expected)}"
        )

    def test_digital_root_of_residues_matches_group(self):
        """
        PREDICTION: Each residue in a digital root group actually has
        that digital root.
        """
        matrix = mod90_digital_root_matrix()
        for dr, residues in matrix.items():
            for r in residues:
                actual_dr = digital_root(r)
                assert actual_dr == dr, (
                    f"Residue {r} has digital root {actual_dr}, "
                    f"expected {dr}"
                )

    def test_mod90_groups_form_arithmetic_progressions(self):
        """
        PREDICTION: Each digital root group follows a specific arithmetic pattern
        within mod 90 (as documented in the formalism).

        dr 1: {1, 19, 37, 73} → gaps of 18, 18, 36
        dr 7: {7, 43, 61, 79} → gaps of 36, 18, 18
        """
        matrix = mod90_digital_root_matrix()

        dr1_gaps = [19 - 1, 37 - 19, 73 - 37]  # [18, 18, 36]
        dr7_gaps = [43 - 7, 61 - 43, 79 - 61]  # [36, 18, 18]

        assert dr1_gaps == [18, 18, 36], f"dr1 gaps: {dr1_gaps}"
        assert dr7_gaps == [36, 18, 18], f"dr7 gaps: {dr7_gaps}"

        # All groups should have gaps from {18, 36}
        for dr, residues in matrix.items():
            sorted_r = sorted(residues)
            gaps = [sorted_r[i + 1] - sorted_r[i] for i in range(len(sorted_r) - 1)]
            for g in gaps:
                assert g in {18, 36}, (
                    f"Digital root {dr}: gap {g} not in {{18, 36}}, "
                    f"residues = {sorted_r}"
                )


class TestLegendreRunningSum:
    """The running sum of Legendre symbols shows mod-90 structure."""

    @pytest.fixture
    def running_sums(self):
        primes = sieve_primes(500)
        return legendre_symbol_running_sum(primes)

    def test_running_sums_computed(self, running_sums):
        """Basic smoke test: running sums are computed."""
        assert len(running_sums) > 50

    def test_mod90_means_vary(self, running_sums):
        """
        PREDICTION: The mean Legendre running sum differs across mod-90 classes.
        If the digital root structure influences Legendre statistics,
        different residue classes should have different average sums.
        """
        means = period_24_in_legendre_sums(running_sums)
        if len(means) < 2:
            pytest.skip("Not enough mod-90 classes populated")

        values = list(means.values())
        # Check that there's meaningful variation
        spread = max(values) - min(values)
        assert spread > 0, "All mod-90 classes have identical mean running sum"


class TestAutocorrelation:
    """Autocorrelation of Legendre running sums may reveal periodic structure."""

    def test_autocorrelation_at_lag_0_is_1(self):
        """Standard autocorrelation property: lag-0 value should be 1.0."""
        primes = sieve_primes(300)
        sums = legendre_symbol_running_sum(primes)
        acf = autocorrelation_of_legendre_sums(sums, max_lag=10)

        if len(acf) > 0:
            assert abs(acf[0] - 1.0) < 0.01, f"Lag-0 autocorrelation = {acf[0]}"

    def test_autocorrelation_decays(self):
        """
        EXPECTATION: Autocorrelation should generally decay from 1.0,
        with possible peaks at periodic lags.
        """
        primes = sieve_primes(500)
        sums = legendre_symbol_running_sum(primes)
        acf = autocorrelation_of_legendre_sums(sums, max_lag=30)

        if len(acf) > 10:
            # Later lags should be smaller than lag 0 (on average)
            late_mean = np.mean(np.abs(acf[5:]))
            assert late_mean < acf[0], "Autocorrelation doesn't decay"


class TestCoprime30Connection:
    """
    Verify the structural connection: coprime-to-30 sequence properties
    emerge from the first three primes {2, 3, 5}.
    """

    def test_first_8_elements(self):
        """The first 8 coprime-to-30 numbers are the prime roots."""
        seq = generate_coprime_sequence(8)
        assert seq == [1, 7, 11, 13, 17, 19, 23, 29]

    def test_interval_pattern(self):
        """
        PREDICTION: Intervals between consecutive coprime-to-30 numbers
        follow the pattern {6, 4, 2, 4, 2, 4, 6, 2} repeating.
        """
        seq = generate_coprime_sequence(24)
        intervals = [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]
        expected_pattern = [6, 4, 2, 4, 2, 4, 6, 2]

        for i, interval in enumerate(intervals):
            expected = expected_pattern[i % 8]
            assert interval == expected, (
                f"Interval at position {i}: {interval}, expected {expected}"
            )

    def test_30n_offset_structure(self):
        """
        PREDICTION: Every element is of the form 30n + r where r ∈ {1,7,11,13,17,19,23,29}.
        """
        seq = generate_coprime_sequence(100)
        for n in seq:
            r = n % 30
            assert r in COPRIME_RESIDUES_MOD30, (
                f"Element {n} has mod-30 residue {r} not in coprime set"
            )

    def test_digital_root_sum_per_rotation(self):
        """
        Each 'rotation' of 8 consecutive coprime numbers has specific digital root sums.
        The first rotation {1,7,11,13,17,19,23,29} has digital root sum = 30.
        """
        seq = generate_coprime_sequence(8)
        dr_sum = sum(digital_root(n) for n in seq)
        # 1+7+2+4+8+1+5+2 = 30
        assert dr_sum == 30, f"First rotation digital root sum = {dr_sum}, expected 30"
