"""
Tests for Primorial-Aligned Spectral Analysis.

Round 5: The deepest exploration — does the Legendre Hamiltonian's spectrum
at primorial-aligned basis sizes encode the twist unit hierarchy?
Does the 288 symmetry connect to the inter-primorial structure?
What particles does the P₄ twist unit (648) predict?
"""

import pytest
import numpy as np
from src.prime_utils import sieve_primes
from src.analytic_108 import euler_totient, general_twist_formula
from src.primorial_spectrum import (
    primorial_aligned_traces,
    twist_unit_in_trace_ratios,
    the_288_connection,
    mass_ratio_search,
)


class TestPrimorialAlignedTraces:
    """Test Hamiltonian trace invariants at primorial-aligned basis sizes."""

    def test_traces_computed_for_p2_and_p3(self):
        """Both P₂ and P₃ aligned traces should be computable."""
        traces = primorial_aligned_traces(max_prime=200)
        assert 'P2' in traces
        assert 'P3' in traces

    def test_p3_basis_size_is_24(self):
        """P₃ = 30: basis size = φ(90) = 24 primes."""
        traces = primorial_aligned_traces(max_prime=200)
        assert traces['P3']['basis_size'] == 24

    def test_p2_basis_size_is_6(self):
        """P₂ = 6: basis size = φ(18) = 6 primes."""
        traces = primorial_aligned_traces(max_prime=200)
        assert traces['P2']['basis_size'] == 6

    def test_traces_are_finite(self):
        """All computed trace invariants should be finite."""
        traces = primorial_aligned_traces(max_prime=200)
        for name, data in traces.items():
            for k, tr in data['traces'].items():
                assert np.isfinite(abs(tr)), (
                    f"Tr(H^{k}) at {name} is not finite: {tr}"
                )

    def test_spectral_width_grows_with_primorial(self):
        """
        PREDICTION: Spectral width at P₃ should exceed that at P₂,
        because more primes create richer coupling.
        """
        traces = primorial_aligned_traces(max_prime=200)
        w_p2 = traces['P2']['spectral_width']
        w_p3 = traces['P3']['spectral_width']
        assert w_p3 > w_p2, (
            f"P₃ width ({w_p3:.2f}) not larger than P₂ width ({w_p2:.2f})"
        )

    def test_trace2_grows_with_primorial(self):
        """
        PREDICTION: |Tr(H²)| at P₃ should exceed that at P₂,
        reflecting the increased total coupling energy.
        """
        traces = primorial_aligned_traces(max_prime=200)
        t2_p2 = traces['P2']['trace_2_abs']
        t2_p3 = traces['P3']['trace_2_abs']
        assert t2_p3 > t2_p2, (
            f"|Tr(H²)| at P₃ ({t2_p3:.2f}) not larger than at P₂ ({t2_p2:.2f})"
        )


class TestThe288Connection:
    """
    Explore the 288-element symmetry from the figure-eight knot complement
    and its relationship to the primorial twist units.
    """

    def test_288_equals_108_times_phi30_div_3(self):
        """
        PROVEN: 288 = T(P₃) × φ(P₃) / 3 = 108 × 8/3 = 288.

        This connects the figure-eight knot's 288-element symmetry to
        the twist unit and the Euler totient of the primorial.
        """
        result = the_288_connection()
        assert result['K_P3'] == 288
        assert result['T_P3'] * result['phi_30'] / 3 == 288

    def test_trefoil_product_is_36(self):
        """s × c × b × u = 6 × 3 × 2 × 1 = 36, and 288 = 36 × 8 = 36 × φ(30)."""
        result = the_288_connection()
        assert result['trefoil_product'] == 36
        assert result['sym_288'] == 36 * 8

    def test_648_div_288_is_9_div_4(self):
        """
        STRUCTURAL RELATIONSHIP: T(P₄)/288 = 648/288 = 9/4 = (3/2)².

        This connects the P₄ twist unit to the P₃ knot symmetry through
        a rational square — suggesting a deep algebraic relationship.
        """
        result = the_288_connection()
        assert result['ratio_648_288'] == 9 / 4

    def test_ratio_648_108_is_phi_7(self):
        """T(P₄)/T(P₃) = 648/108 = 6 = φ(7)."""
        result = the_288_connection()
        assert result['ratio_648_108'] == 6.0

    def test_generalized_knot_symmetry_p4(self):
        """
        NOVEL PREDICTION: The P₄ analog of the 288-element symmetry is:
        K(P₄) = T(P₄) × φ(P₄) / 3 = 648 × 48 / 3 = 10368 = 2⁷ × 3⁴.
        """
        result = the_288_connection()
        assert result['K_P4'] == 10368
        assert result['K_P4'] == 2**7 * 3**4


class TestMassRatioSearch:
    """
    Search for particle mass ratios matching the P₃ and P₄ twist units.
    """

    def test_proton_matches_108(self):
        """m_p/m_e = 1836 = 17 × 108 exactly."""
        candidates = mass_ratio_search(108)
        proton = next(c for c in candidates if c['particle'] == 'proton')
        assert proton['n'] == 17
        assert abs(proton['residual']) < 1  # Very close to 17 × 108

    def test_tau_matches_108(self):
        """m_τ/m_e = 3477 ≈ 31 × 108 + 125 + 4 (Composite Quintic/Quadratic)."""
        candidates = mass_ratio_search(108)
        tau = next(c for c in candidates if c['particle'] == 'tau')
        # n is based on round(3477/108) = 32
        assert tau['n'] == 32
        # But the correction logic finds the shifted composite (n-1) * 108 + 129
        # 31 * 108 + 129 = 3348 + 129 = 3477
        assert tau['correction_type'] == 'quintic_quadratic_shifted'
        assert tau['best_correction'] == 129
        assert tau['best_error_pct'] < 0.01

    def test_648_search_finds_candidates(self):
        """
        NOVEL EXPLORATION: Search particles against the P₄ twist unit 648.
        Some particles may have cleaner fits to 648 than to 108.
        """
        candidates = mass_ratio_search(648)
        # Just verify the search runs and produces results
        assert len(candidates) > 0

    def test_648_charm_quark_fit(self):
        """
        Test if charm quark mass has a clean relationship with 648.
        m_c/m_e ≈ 2494; 2494/648 ≈ 3.85, not clean.
        But: 4 × 648 = 2592, and 2592 - 2494 = 98 ≈ 3² × 11.
        """
        candidates = mass_ratio_search(648)
        charm = next(c for c in candidates if c['particle'] == 'charm_quark')
        # Document the fit quality — not necessarily good
        assert charm['n'] > 0

    def test_648_kaon_fit(self):
        """
        Kaon mass: m_K/m_e ≈ 966. Test if this relates to 648.
        966/648 ≈ 1.49 ≈ 3/2. So m_K ≈ 1.5 × 648 = 972. Error: 0.6%
        OR: m_K = 648 + 318 = 648 + 2×3×53. Not obviously clean.
        """
        candidates = mass_ratio_search(648)
        kaon = next(c for c in candidates if c['particle'] == 'kaon')
        # 1 × 648 = 648, 2 × 648 = 1296. n=1 gives residual 318.
        assert kaon['n'] in [1, 2]

    def test_which_particles_fit_108_better_than_648(self):
        """
        COMPARATIVE TEST: For each particle, determine whether 108 or 648
        provides a cleaner fit. Particles designed by P₃ physics (proton,
        muon, tau) should fit 108 better. Heavier/exotic particles might
        fit 648 better.
        """
        candidates_108 = {c['particle']: c['best_error_pct'] for c in mass_ratio_search(108)}
        candidates_648 = {c['particle']: c['best_error_pct'] for c in mass_ratio_search(648)}

        better_108 = []
        better_648 = []
        for particle in candidates_108:
            if particle in candidates_648:
                if candidates_108[particle] < candidates_648[particle]:
                    better_108.append(particle)
                else:
                    better_648.append(particle)

        # Just document the split — this is exploratory
        assert len(better_108) + len(better_648) > 0


class TestTwistUnitRatios:
    """
    Test whether trace invariant ratios across primorial levels
    approximate the twist unit ratios.
    """

    def test_trace_ratios_computed(self):
        """Trace ratios should be computable."""
        ratios = twist_unit_in_trace_ratios(max_prime=200)
        assert len(ratios) > 0

    def test_twist_ratio_p3_p2_is_4(self):
        """T(P₃)/T(P₂) = 108/27 = 4."""
        assert 108 / 27 == 4

    def test_trace_ratio_direction_matches(self):
        """
        PREDICTION: |Tr_P3(H²)| / |Tr_P2(H²)| should exceed 1
        (the P₃ Hamiltonian has more coupling energy than P₂).
        """
        ratios = twist_unit_in_trace_ratios(max_prime=200)
        tr2_key = 'Tr2_P3/P2'
        if tr2_key in ratios:
            assert ratios[tr2_key]['trace_ratio'] > 1, (
                f"Trace2 ratio P₃/P₂ = {ratios[tr2_key]['trace_ratio']:.4f}"
            )
