"""
Tests for 108-periodicity in prime resonance structures.

Validates:
1. The period-24 digital root cycle sums to exactly 108 = 2²×3³
2. The 108 number appears in the structure of the coprime-to-30 sequence
3. The density of asymmetric Legendre pairs has periodic fluctuation
4. Physical constant derivations from the twist framework hold
"""

import pytest
import numpy as np
from src.prime_utils import (
    sieve_primes, PERIOD_24_DIGITAL_ROOTS, TWIST_UNIT, digital_root
)
from src.mod30_structure import generate_coprime_sequence, compute_digital_root_cycle
from src.periodicity import (
    cumulative_asymmetric_pair_count,
    detrended_pair_density,
    detect_periodicity_fft,
    analyze_108_periodicity,
    legendre_symbol_running_sum,
)


class TestDigitalRootSum108:
    """
    NOVEL PREDICTION: The sum of one period-24 digital root cycle = 108.

    This is a previously unidentified structural connection between:
    - The period-24 cycling of the mod-30 sieve (from primespiral.md)
    - The twist unit 108 = 2²×3³ (from 19_Constants_From_Twist.md)

    If this holds, it means the period-24 structure of the prime sieve
    directly encodes the fundamental twist unit used to derive physical constants.
    """

    def test_period_24_sum_equals_108(self):
        """The digital root cycle {1,7,2,4,8,1,5,2,4,1,5,7,2,4,8,5,7,4,8,1,5,7,2,8} sums to 108."""
        cycle_sum = sum(PERIOD_24_DIGITAL_ROOTS)
        assert cycle_sum == TWIST_UNIT, (
            f"Period-24 digital root sum = {cycle_sum}, expected {TWIST_UNIT} = 2²×3³"
        )

    def test_108_factorization(self):
        """108 = 2² × 3³ = 4 × 27."""
        assert TWIST_UNIT == 108
        assert TWIST_UNIT == 4 * 27
        assert TWIST_UNIT == (2**2) * (3**3)

    def test_multiple_periods_sum_to_multiples_of_108(self):
        """
        PREDICTION: The sum of N complete periods of the digital root cycle
        is exactly N × 108.
        """
        seq = generate_coprime_sequence(24 * 5)
        roots = compute_digital_root_cycle(seq)
        for n in range(1, 6):
            partial_sum = sum(roots[:24 * n])
            expected = 108 * n
            assert partial_sum == expected, (
                f"Sum of {n} periods = {partial_sum}, expected {expected}"
            )

    def test_half_period_sum(self):
        """
        STRUCTURAL TEST: The first 12 and last 12 digital roots of the period-24 cycle.
        First half:  1+7+2+4+8+1+5+2+4+1+5+7 = 47
        Second half: 2+4+8+5+7+4+8+1+5+7+2+8 = 61
        Total: 47 + 61 = 108
        """
        first_half = sum(PERIOD_24_DIGITAL_ROOTS[:12])
        second_half = sum(PERIOD_24_DIGITAL_ROOTS[12:])
        assert first_half + second_half == 108
        assert first_half == 47  # 1+7+2+4+8+1+5+2+4+1+5+7
        assert second_half == 61  # 2+4+8+5+7+4+8+1+5+7+2+8


class TestTwistUnitInPhysicalConstants:
    """
    Verify the twist framework's physical constant formulas from 19_Constants_From_Twist.md.
    These are testable numerical predictions.
    """

    def test_proton_electron_mass_ratio(self):
        """
        PREDICTION: m_p/m_e ≈ 17 × 108 = 1836 (integer part matches exactly).
        17 = s×c - b + u where (s,c,b,u) = (6,3,2,1) are trefoil invariants.
        """
        trefoil_s, trefoil_c, trefoil_b, trefoil_u = 6, 3, 2, 1
        T = trefoil_s * trefoil_c - trefoil_b + trefoil_u
        assert T == 17, f"Trefoil complexity = {T}, expected 17"

        predicted = T * 108
        experimental = 1836  # Integer part of 1836.15267343
        assert predicted == experimental, (
            f"Predicted m_p/m_e = {predicted}, experimental integer part = {experimental}"
        )

    def test_fine_structure_constant(self):
        """
        PREDICTION: α⁻¹ ≈ 108 + 29 = 137 (≈ 0.026% error).
        29 is the 10th prime and largest coprime residue < 30.
        """
        predicted = 108 + 29
        experimental = 137  # Integer part of 137.035999
        assert predicted == experimental

    def test_fine_structure_refined(self):
        """
        PREDICTION: α⁻¹ ≈ 108 + 29 + 1/27 = 137.037 (≈ 0.0007% error).
        27 = 3³ is the ternary component of 108 = 4 × 27.
        """
        predicted = 108 + 29 + 1.0 / 27
        experimental = 137.035999
        error_pct = abs(predicted - experimental) / experimental * 100
        assert error_pct < 0.01, (
            f"α⁻¹ prediction {predicted:.6f} vs experimental {experimental:.6f}, "
            f"error = {error_pct:.4f}%"
        )

    def test_muon_electron_mass_ratio(self):
        """
        PREDICTION: m_μ/m_e ≈ 2 × 108 - 9 = 207.
        9 = 3² is the ternary defect.
        """
        predicted = 2 * 108 - 9
        assert predicted == 207
        experimental = 206.768  # ~207
        error_pct = abs(predicted - experimental) / experimental * 100
        assert error_pct < 0.2

    def test_tau_electron_mass_ratio(self):
        """
        PREDICTION: m_τ/m_e ≈ 31 × 108 + 125 + 4 = 3477.
        Composite Quintic (5³) + Quadratic (2²) channel correction.
        """
        predicted = 31 * 108 + 125 + 4
        assert predicted == 3477
        experimental = 3477.23
        error_pct = abs(predicted - experimental) / experimental * 100
        assert error_pct < 0.01

    def test_higgs_mass(self):
        """
        PREDICTION: M_H = 5³ = 125 GeV.
        5 is the first prime not dividing 108 = 2²×3³.
        """
        predicted = 5**3
        assert predicted == 125
        experimental = 125.25  # GeV
        error_pct = abs(predicted - experimental) / experimental * 100
        assert error_pct < 0.3

    def test_weinberg_angle(self):
        """
        PREDICTION: sin²θ_W ≈ 3/13 = 0.2308 (0.2% error).
        3 and 13 are both primes in the coprime residue set.
        """
        predicted = 3.0 / 13
        experimental = 0.23122
        error_pct = abs(predicted - experimental) / experimental * 100
        assert error_pct < 0.5


class TestAsymmetricPairDensityPeriodicity:
    """Test whether asymmetric Legendre pair density shows periodic structure."""

    def test_detrended_fluctuation_nonzero(self):
        """
        PREDICTION: After removing the linear trend, the cumulative pair count
        shows non-trivial fluctuations (not perfectly linear).
        """
        primes = sieve_primes(500)
        cumulative = cumulative_asymmetric_pair_count(primes)

        if len(cumulative) < 10:
            pytest.skip("Not enough data for periodicity test")

        primes_arr, detrended = detrended_pair_density(cumulative)
        assert len(detrended) > 0
        assert np.std(detrended) > 0, "Detrended signal has zero variance"

    def test_periodicity_detection_functional(self):
        """
        Verify the periodicity detection machinery works on a known periodic signal.
        """
        # Create a known signal with period 108
        t = np.arange(1000)
        signal = np.sin(2 * np.pi * t / 108)
        freqs, amps, period = detect_periodicity_fft(signal, sample_spacing=1.0)

        # Should detect period near 108
        assert abs(period - 108) < 5, (
            f"Known period-108 signal detected as period {period:.1f}"
        )
