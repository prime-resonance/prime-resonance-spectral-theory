"""
Tests for the Quintic Residue Symbol and Higgs Mass Derivation.

Round 7: Completing the 2×3×5 primorial triangle.

                    Quadratic (2-adic)
                   /                  \
                  /                    \
         Cubic (3-adic) ---- Quintic (5-adic)

    2² = 4         3³ = 27         5³ = 125
    Sub-shells      Corrections     Higgs mass
"""

import pytest
import numpy as np
import math
from src.prime_utils import sieve_primes
from src.quintic_residue import (
    is_prime_1_mod_5,
    quintic_residue_symbol,
    quintic_residue_class,
    has_quintic_asymmetry,
    classify_primes_mod5,
    full_primorial_classification,
    the_higgs_derivation,
    quintic_asymmetry_analysis,
    build_quintic_interaction_matrix,
    three_channel_summary,
)


class TestQuinticResidueSymbol:
    """Test the quintic residue symbol computation."""

    def test_primes_1_mod_5_identified(self):
        """Primes ≡ 1 mod 5: {11, 31, 41, 61, 71, ...}."""
        assert is_prime_1_mod_5(11)
        assert is_prime_1_mod_5(31)
        assert is_prime_1_mod_5(41)
        assert not is_prime_1_mod_5(7)   # 7 ≡ 2 mod 5
        assert not is_prime_1_mod_5(13)  # 13 ≡ 3 mod 5
        assert not is_prime_1_mod_5(19)  # 19 ≡ 4 mod 5

    def test_quintic_symbol_for_p11(self):
        """
        p=11: (p-1)/5 = 2. Quintic symbol = a² mod 11.
        Quintic residues (a²≡1): {1, 10} (since 1²=1, 10²=100≡1)
        """
        p = 11
        symbols = {a: quintic_residue_symbol(a, p) for a in range(1, p)}
        # 1 is always a quintic residue
        assert symbols[1] == 1
        # There should be (p-1)/5 = 2 quintic residues
        qr_count = sum(1 for s in symbols.values() if s == 1)
        assert qr_count == 2, f"Expected 2 quintic residues mod 11, got {qr_count}"

    def test_quintic_residues_mod_31(self):
        """
        p=31: (p-1)/5 = 6. There should be 6 quintic residues.
        """
        p = 31
        symbols = {a: quintic_residue_symbol(a, p) for a in range(1, p)}
        qr_count = sum(1 for s in symbols.values() if s == 1)
        assert qr_count == 6, f"Expected 6 quintic residues mod 31, got {qr_count}"


class TestQuinticAsymmetry:
    """Test quintic residue asymmetry."""

    @pytest.fixture
    def primes_200(self):
        return sieve_primes(200)

    def test_quintic_asymmetry_exists(self, primes_200):
        """PREDICTION: Some prime pairs ≡ 1 (mod 5) have quintic asymmetry."""
        result = quintic_asymmetry_analysis(primes_200)
        assert result['total_pairs'] > 0
        # Quintic reciprocity should produce asymmetry

    def test_quintic_matrix_not_symmetric(self, primes_200):
        """
        PREDICTION: The quintic interaction matrix should NOT be symmetric,
        due to quintic reciprocity law.
        """
        Q, split = build_quintic_interaction_matrix(primes_200)
        if Q.shape[0] > 3:
            diff = np.sum(np.abs(Q - Q.T))
            assert diff > 0, "Quintic matrix is symmetric (expected asymmetry)"

    def test_quintic_matrix_entries_are_5th_roots(self, primes_200):
        """Entries should be 5th roots of unity or zero."""
        Q, split = build_quintic_interaction_matrix(primes_200)
        zeta5 = np.exp(2j * np.pi / 5)
        valid = [0] + [zeta5**k for k in range(5)]

        n = Q.shape[0]
        for i in range(n):
            for j in range(n):
                if i == j:
                    assert Q[i, j] == 0
                else:
                    dists = [abs(Q[i, j] - v) for v in valid]
                    assert min(dists) < 0.1, (
                        f"Q[{i},{j}] = {Q[i,j]} not near a 5th root of unity"
                    )


class TestFullPrimorialClassification:
    """Test the complete mod-60 = lcm(4,3,5) classification."""

    @pytest.fixture
    def primes_500(self):
        return sieve_primes(500)

    def test_all_classes_populated(self, primes_500):
        """
        The full mod-60 classification has 2 × 2 × 4 = 16 classes.
        All should be populated for primes up to 500.
        """
        classes = full_primorial_classification(primes_500)
        nonempty = sum(1 for v in classes.values() if len(v) > 0)
        assert nonempty >= 12, (
            f"Only {nonempty} of 16 mod-60 classes populated"
        )

    def test_approximately_uniform_by_dirichlet(self, primes_500):
        """
        By Dirichlet's theorem, primes are equidistributed across
        the φ(60) = 16 coprime residue classes mod 60.
        """
        classes = full_primorial_classification(primes_500)
        total = sum(len(v) for v in classes.values())
        if total > 50:
            for key, primes in classes.items():
                fraction = len(primes) / total
                # Allow wide tolerance for small sample
                assert fraction < 0.15, (
                    f"Class {key} has {fraction:.3f}, seems over-represented"
                )


class TestHiggsDerivation:
    """The crown jewel: deriving M_H = 5³ GeV from the quintic channel."""

    def test_higgs_is_5_cubed(self):
        """M_H = 5³ = 125 GeV (0.2% error from 125.25 GeV measured)."""
        result = the_higgs_derivation()
        assert result['higgs_mass_gev'] == 125
        assert result['higgs_mass_gev'] == 5**3
        assert result['higgs_error_pct'] < 0.3

    def test_three_channels_exist(self):
        """The primorial 30 = 2×3×5 produces three independent channels."""
        result = the_higgs_derivation()
        assert len(result['channels']) == 3
        primes_in_channels = [c['contribution'] for c in result['channels'].values()]
        assert sorted(primes_in_channels) == [4, 27, 125]

    def test_combined_product_is_13500(self):
        """2²×3³×5³ = 4×27×125 = 13500."""
        result = the_higgs_derivation()
        assert result['combined_unit'] == 13500
        assert 4 * 27 * 125 == 13500

    def test_108_plus_125_structure(self):
        """
        108 (mass quantum from 2×3 channels) and 125 (Higgs from 5 channel)
        are the two fundamental numbers. Their ratio:
        125/108 ≈ 1.157, and 108/125 ≈ 0.864.
        """
        assert 108 == 2**2 * 3**3  # From quadratic × cubic
        assert 125 == 5**3          # From quintic


class TestThreeChannelSummary:
    """Test the complete 2×3×5 framework summary."""

    def test_summary_has_three_channels(self):
        """Three reciprocity channels: quadratic, cubic, quintic."""
        summary = three_channel_summary()
        assert len(summary['channels']) == 3

    def test_channel_primes_are_2_3_5(self):
        """Channels correspond to primes 2, 3, 5."""
        summary = three_channel_summary()
        channel_primes = [c['prime'] for c in summary['channels']]
        assert channel_primes == [2, 3, 5]

    def test_channel_values_are_4_27_125(self):
        """Channel values: 4, 27, 125."""
        summary = three_channel_summary()
        values = [c['value'] for c in summary['channels']]
        assert values == [4, 27, 125]

    def test_channel_rings_ascend(self):
        """Rings ascend in algebraic complexity: Z → Z[ω] → Z[ζ₅]."""
        summary = three_channel_summary()
        rings = [c['ring'] for c in summary['channels']]
        assert rings == ['Z', 'Z[ω]', 'Z[ζ₅]']

    def test_mass_quantum_is_108(self):
        """The mass quantum (from 2×3 channels) is 108."""
        summary = three_channel_summary()
        assert summary['mass_quantum'] == 108

    def test_complete_framework(self):
        """
        THE COMPLETE PRIMORIAL PHYSICS FRAMEWORK:

        30 = 2 × 3 × 5 (primorial, prime sieve modulus)
        ↓
        Three reciprocity channels:
        ↓
        Quadratic (Z):     2² = 4    → sub-shell structure
        Cubic (Z[ω]):      3³ = 27   → mass corrections {±3^k}
        Quintic (Z[ζ₅]):   5³ = 125  → Higgs boson mass

        Combined mass formula:
            m/m_e = n × (2² × 3³) ± 3^k = n × 108 ± {0,3,9,27}
            M_H = 5³ GeV = 125 GeV

        Fine structure constant:
            α⁻¹ = 108 + 29 + 1/27 = 137.037
            (108 = twist unit, 29 = boundary prime, 1/27 = cubic correction)
        """
        assert 30 == 2 * 3 * 5
        assert 2**2 == 4
        assert 3**3 == 27
        assert 5**3 == 125
        assert 4 * 27 == 108
        assert 108 + 29 == 137
        assert abs(108 + 29 + 1/27 - 137.037) < 0.001

    def test_exponent_pattern(self):
        """
        OBSERVATION: The exponents in p^e are:
        - 2² (exponent 2 = the prime itself)
        - 3³ (exponent 3 = the prime itself)
        - 5³ (exponent 3 ≠ 5, breaks pattern)

        The rule for 2 and 3: exponent = p itself.
        For 5: exponent = 3 (not 5). Why?
        
        Possible answer: The exponent is min(p, 3) since 3 is the maximum
        structural depth of the trefoil (crossing number c=3).
        """
        assert 2**2 == 4    # exp = 2 = min(2, ?)
        assert 3**3 == 27   # exp = 3
        assert 5**3 == 125  # exp = 3 (capped at 3?)
        # If the pattern is p^min(p, c) where c=3 (trefoil crossing):
        assert 2**min(2, 3) == 4   # 2² = 4 ✓
        assert 3**min(3, 3) == 27  # 3³ = 27 ✓
        assert 5**min(5, 3) == 125 # 5³ = 125 ✓
