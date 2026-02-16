"""
Tests for Mixing Angles and P4 Scale Physics.
"""

import pytest
import math
from src.mixing_angles import mixing_angle_search, p4_scale_physics

class TestMixingAngles:
    """
    Test predictions for CKM and PMNS mixing angles.
    """
    
    @pytest.fixture
    def angles(self):
        return mixing_angle_search()
    
    def test_weinberg_angle_prediction(self, angles):
        """
        PREDICTION: sin²θ_W ≈ 3/13.
        """
        assert angles['weinberg']['error_pct'] < 0.5, \
            f"Weinberg angle error {angles['weinberg']['error_pct']:.2f}% too high"

    def test_cabibbo_angle_prediction(self, angles):
        """
        PREDICTION: sin(θ_C) ≈ 29/128.
        """
        assert angles['cabibbo']['error_pct'] < 1.0, \
            f"Cabibbo angle error {angles['cabibbo']['error_pct']:.2f}% too high"
            
    def test_pmns_12_prediction(self, angles):
        """
        PREDICTION: sin²θ_12 ≈ 4/13.
        """
        assert angles['pmns_12']['error_pct'] < 2.0, \
            f"PMNS θ_12 error {angles['pmns_12']['error_pct']:.2f}% too high"

    def test_pmns_23_prediction(self, angles):
        """
        PREDICTION: sin²θ_23 ≈ 4/7.
        """
        assert angles['pmns_23']['error_pct'] < 5.0, \
            f"PMNS θ_23 error {angles['pmns_23']['error_pct']:.2f}% too high"
            
    def test_pmns_13_prediction(self, angles):
        """
        PREDICTION: sin²θ_13 ≈ 1/45.
        """
        assert angles['pmns_13']['error_pct'] < 5.0, \
            f"PMNS θ_13 error {angles['pmns_13']['error_pct']:.2f}% too high"

class TestP4ScalePhysics:
    """
    Test predictions for P4 scale physics.
    """
    
    @pytest.fixture
    def scales(self):
        return p4_scale_physics()
        
    def test_k_p4_mass_b_meson(self, scales):
        """
        PREDICTION: K(P4) * m_e ≈ 5.3 GeV, corresponding to B meson mass range.
        B mesons are ~5279 MeV.
        """
        predicted = scales['K_P4_mass']
        # 10368 * 0.511 = 5298.048 MeV
        experimental = 5279.0
        error = abs(predicted - experimental) / experimental * 100
        assert error < 1.0, f"B meson mass prediction error {error:.2f}% too high"

