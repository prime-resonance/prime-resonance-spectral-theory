"""
Tests for Neutrino Mass Predictions from the Inverse Twist Tower.

Round 8: The framework's most falsifiable prediction — neutrino masses
derived from the INVERSE of the twist unit hierarchy.

If positive masses are n × 108 ± 3^k,
then neutrino masses are 1/(108² × 3^k) for k ∈ {0, 1, 2}.
"""

import pytest
import math
from src.neutrino_prediction import (
    inverse_twist_mass_ratio,
    predict_neutrino_masses,
    mass_hierarchy_table,
    neutrino_mixing_from_3adic,
    M_ELECTRON_EV,
)


class TestInverseTwistTower:
    """Test the inverse twist mass formula."""

    def test_inverse_level_0(self):
        """1/108⁰ = 1 (electron mass itself)."""
        assert inverse_twist_mass_ratio(0, 0) == 1.0

    def test_inverse_level_1(self):
        """1/108 ≈ 9.26 × 10⁻³."""
        ratio = inverse_twist_mass_ratio(1, 0)
        assert abs(ratio - 1/108) < 1e-10

    def test_inverse_level_2(self):
        """1/108² ≈ 8.57 × 10⁻⁵."""
        ratio = inverse_twist_mass_ratio(2, 0)
        assert abs(ratio - 1/11664) < 1e-15

    def test_3adic_correction_k1(self):
        """1/(108² × 3) ≈ 2.86 × 10⁻⁵."""
        ratio = inverse_twist_mass_ratio(2, 1)
        assert abs(ratio - 1/34992) < 1e-15

    def test_3adic_correction_k2(self):
        """1/(108² × 9) ≈ 9.53 × 10⁻⁶."""
        ratio = inverse_twist_mass_ratio(2, 2)
        assert abs(ratio - 1/104976) < 1e-15


class TestNeutrinoPredictions:
    """Test the concrete neutrino mass predictions."""

    @pytest.fixture
    def predictions(self):
        return predict_neutrino_masses()

    def test_three_mass_eigenstates(self, predictions):
        """Three neutrino mass eigenstates predicted."""
        assert 'ν_1' in predictions
        assert 'ν_2' in predictions
        assert 'ν_3' in predictions

    def test_normal_ordering(self, predictions):
        """
        PREDICTION: Normal ordering m₁ < m₂ < m₃.
        This comes from the 3-adic tower: k=2 → lightest, k=0 → heaviest.
        """
        m1 = predictions['ν_1']['mass_eV']
        m2 = predictions['ν_2']['mass_eV']
        m3 = predictions['ν_3']['mass_eV']
        assert m1 < m2 < m3, f"Not normal ordering: {m1:.3e} < {m2:.3e} < {m3:.3e}"

    def test_mass_sum_below_cosmological_bound(self, predictions):
        """
        PREDICTION: Σm_ν = 73.2 meV = 0.073 eV < 0.12 eV (Planck 2018 bound).
        Formula: m_ν = m_e / (108³ × 8 × 3^k) for k ∈ {0,1,2}.
        """
        total = predictions['mass_sum_eV']
        bound = predictions['experimental']['sum_bound']
        assert total < bound, (
            f"Predicted Σm_ν = {total:.4f} eV exceeds Planck bound {bound} eV"
        )

    def test_heaviest_mass_near_50_meV(self, predictions):
        """
        PREDICTION: m₃ = m_e / (108³ × 8) ≈ 50.7 meV ≈ 0.0507 eV.
        This matches √Δm²₃₂ ≈ 49.5 meV at 2.4% error.
        """
        m3 = predictions['ν_3']['mass_eV']
        assert 0.01 < m3 < 0.1, f"m₃ = {m3:.5f} eV, expected 0.01-0.1 eV"
        # Compare to √Δm²₃₂ ≈ 0.050 eV
        error_pct = abs(m3 - 0.0495) / 0.0495 * 100
        assert error_pct < 10, f"m₃ = {m3:.5f} eV, {error_pct:.1f}% from ~49.5 meV"

    def test_atmospheric_mass_squared_difference(self, predictions):
        """
        PREDICTION: Δm²₃₂ = m₃² - m₂² ≈ 2.285 × 10⁻³ eV².
        Measured: 2.453 × 10⁻³ eV². Error: 6.8%.

        THIS IS THE STRONGEST NEUTRINO PREDICTION: the atmospheric 
        mass-squared difference is predicted at 6.8% accuracy from
        pure number theory (108³ × φ(30) × 3^k).
        """
        dm32_predicted = predictions['dm32_squared']
        dm32_measured = predictions['experimental']['dm32_sq_measured']
        
        error_pct = abs(dm32_predicted - dm32_measured) / dm32_measured * 100
        assert error_pct < 10, (
            f"Δm²₃₂: predicted {dm32_predicted:.3e}, measured {dm32_measured:.3e}, "
            f"error {error_pct:.1f}%"
        )

    def test_mass_squared_difference_solar(self, predictions):
        """
        PREDICTION: Δm²₂₁ should be smaller than Δm²₃₂.
        """
        dm21 = predictions['dm21_squared']
        dm32 = predictions['dm32_squared']
        assert dm21 < dm32, (
            f"Δm²₂₁ ({dm21:.3e}) should be < Δm²₃₂ ({dm32:.3e})"
        )


class TestMassHierarchyCompleteness:
    """Test the complete Standard Model mass hierarchy."""

    def test_hierarchy_spans_12_orders(self):
        """The mass table should span ~12 orders of magnitude."""
        table = mass_hierarchy_table()
        lightest = table[0]['ratio']
        heaviest = table[-1]['ratio']
        span = math.log10(heaviest / lightest)
        assert span > 10, f"Hierarchy spans only {span:.1f} orders"

    def test_all_particles_present(self):
        """At least 15 particles in the hierarchy."""
        table = mass_hierarchy_table()
        assert len(table) >= 15

    def test_neutrinos_are_lightest(self):
        """Neutrinos should be at the bottom of the hierarchy."""
        table = mass_hierarchy_table()
        assert table[0]['name'].startswith('ν')
        assert table[1]['name'].startswith('ν')
        assert table[2]['name'].startswith('ν')

    def test_top_quark_is_heaviest(self):
        """Top quark should be near the top of the hierarchy."""
        table = mass_hierarchy_table()
        assert table[-1]['name'] == 't'


class TestNeutrinoMixing:
    """Test predictions from the 3-adic mixing structure."""

    def test_consecutive_ratio_is_3(self):
        """
        FIRM PREDICTION: m₃/m₂ = 3 and m₂/m₁ = 3, exactly.
        This is a direct consequence of the 3-adic correction tower.
        """
        mixing = neutrino_mixing_from_3adic()
        assert mixing['m3_to_m2_ratio'] == 3.0
        assert mixing['m2_to_m1_ratio'] == 3.0

    def test_mass_squared_ratio(self):
        """
        DERIVED PREDICTION: The ratio Δm²₂₁ / Δm²₃₂.
        
        With m₁:m₂:m₃ = 1:3:9 (geometric with ratio 3):
        Δm²₂₁ = m₂² - m₁² = 9m₁² - m₁² = 8m₁²
        Δm²₃₂ = m₃² - m₂² = 81m₁² - 9m₁² = 72m₁²
        Ratio: 8/72 = 1/9

        Experimental: Δm²₂₁/Δm²₃₂ ≈ 7.53×10⁻⁵/2.453×10⁻³ ≈ 0.031 = ~1/33
        
        Our prediction: 1/9 ≈ 0.111.
        This is off by factor ~3.6 from experiment.
        """
        mixing = neutrino_mixing_from_3adic()
        predicted_ratio = mixing['dm21_sq_over_dm32_sq']
        
        # The prediction is 1/9 = 0.111
        assert abs(predicted_ratio - 1/9) < 0.001, (
            f"Mass-squared difference ratio = {predicted_ratio:.4f}, expected 1/9 = 0.111"
        )
        
        # For comparison with experiment
        experimental_ratio = 7.53e-5 / 2.453e-3  # ≈ 0.031
        # Our prediction (0.111) is ~3.6× too large, suggesting the actual
        # spacing might be a DIFFERENT power of 3, not consecutive

    def test_geometric_spacing(self):
        """The mass eigenvalues are geometrically spaced with ratio 3."""
        mixing = neutrino_mixing_from_3adic()
        assert mixing['geometric_spacing'] == 3.0
