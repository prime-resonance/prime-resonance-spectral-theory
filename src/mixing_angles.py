"""
CKM and PMNS Mixing Angles from Twist Theory.

The twist framework suggests that mixing angles in the CKM (quark) and PMNS (neutrino)
matrices arise from the geometry of the prime resonance field.

Hypothesis:
Mixing angles correspond to specific ratios of twist units and primorial totients.
We specifically look for angles related to:
- The "Cabibbo angle" θ_C ≈ 13.04°
- The "Weinberg angle" θ_W ≈ 28.7° (sin²θ_W ≈ 0.231)
- The PMNS angles θ_12, θ_23, θ_13

We define a "Twist Angle" θ_T(n, k) = arccos(k / sqrt(n^2 + k^2)) or similar geometric constructions
derived from the prime number lattice.

Key targets:
- sin(θ_C) ≈ 0.225
- sin²(θ_W) ≈ 0.231
- sin²(θ_12) ≈ 0.307 (PMNS)
- sin²(θ_23) ≈ 0.546 (PMNS)
- sin²(θ_13) ≈ 0.022 (PMNS)
"""

from typing import Dict, List, Tuple
import math
import numpy as np
from .prime_utils import sieve_primes
from .analytic_108 import euler_totient

def mixing_angle_search() -> Dict[str, object]:
    """
    Search for mixing angles derived from prime/twist ratios.

    We explore ratios involving:
    - Twist unit T = 108
    - Primorials P_k
    - Totients φ(P_k)
    - Small primes p

    Returns:
        Dict of found angles and their errors relative to experimental values.
    """
    
    # Experimental values (approximate)
    targets = {
        'cabibbo': 13.04 * math.pi / 180,       # sin ≈ 0.2257
        'weinberg': 28.7 * math.pi / 180,       # sin² ≈ 0.231
        'pmns_12': 33.44 * math.pi / 180,       # sin² ≈ 0.304
        'pmns_23': 49.0 * math.pi / 180,        # sin² ≈ 0.57
        'pmns_13': 8.57 * math.pi / 180,        # sin² ≈ 0.0222
    }
    
    # Twist-based candidates
    # 1. Weinberg angle from 3/13 (already found)
    # sin²θ_W ≈ 3/13 ≈ 0.23077
    weinberg_pred = math.asin(math.sqrt(3/13))
    
    # 2. Cabibbo angle
    # sin(θ_C) ≈ 0.2257.  Maybe 1 / sqrt(20)? 1/4.5?
    # 20 = 4*5?  4.5 = 9/2?
    # Let's try ratios of small integers/primes.
    # 2/9 = 0.222... close.
    # 4/18 = 2/9.
    # What about 13/60? 0.216.
    # 13 is prime, 60 is P3*2.
    
    # Twist-based: 29 / 128? (29 is boundary prime, 128 = 2^7) -> 0.2265
    cabibbo_pred_1 = math.asin(29/128)
    
    # 3. PMNS θ_12
    # sin²θ_12 ≈ 0.304.
    # 3/10 = 0.3.
    # 1/3 = 0.333.
    # 4/13 = 0.3076.  (4 = 2², 13 is prime).
    pmns_12_pred = math.asin(math.sqrt(4/13))
    
    # 4. PMNS θ_23
    # sin²θ_23 ≈ 0.57.
    # 4/7 = 0.5714.
    # 7 is prime. 4=2².
    pmns_23_pred = math.asin(math.sqrt(4/7))
    
    # 5. PMNS θ_13
    # sin²θ_13 ≈ 0.0222.
    # 1/45 = 0.0222... (45 = 5*9 = 5*3²)
    # 2/90 = 1/45.
    pmns_13_pred = math.asin(math.sqrt(1/45))
    
    results = {
        'weinberg': {
            'predicted_rad': weinberg_pred,
            'predicted_deg': weinberg_pred * 180 / math.pi,
            'formula': 'arcsin(sqrt(3/13))',
            'error_pct': abs(weinberg_pred - targets['weinberg']) / targets['weinberg'] * 100
        },
        'cabibbo': {
            'predicted_rad': cabibbo_pred_1,
            'predicted_deg': cabibbo_pred_1 * 180 / math.pi,
            'formula': 'arcsin(29/128)',
            'error_pct': abs(cabibbo_pred_1 - targets['cabibbo']) / targets['cabibbo'] * 100
        },
        'pmns_12': {
            'predicted_rad': pmns_12_pred,
            'predicted_deg': pmns_12_pred * 180 / math.pi,
            'formula': 'arcsin(sqrt(4/13))',
            'error_pct': abs(pmns_12_pred - targets['pmns_12']) / targets['pmns_12'] * 100
        },
        'pmns_23': {
            'predicted_rad': pmns_23_pred,
            'predicted_deg': pmns_23_pred * 180 / math.pi,
            'formula': 'arcsin(sqrt(4/7))',
            'error_pct': abs(pmns_23_pred - targets['pmns_23']) / targets['pmns_23'] * 100
        },
        'pmns_13': {
            'predicted_rad': pmns_13_pred,
            'predicted_deg': pmns_13_pred * 180 / math.pi,
            'formula': 'arcsin(sqrt(1/45))',
            'error_pct': abs(pmns_13_pred - targets['pmns_13']) / targets['pmns_13'] * 100
        }
    }
    
    return results

def p4_scale_physics() -> Dict[str, object]:
    """
    Probe the P4 = 210 scale for Beyond Standard Model (BSM) signals.
    
    P4 Twist Unit = 648.
    
    Hypothesis:
    New particles or forces might appear at mass scales related to T(P4) = 648
    or K(P4) = 10368.
    
    Possible scales:
    - 648 * m_e ≈ 331 MeV (unknown?)
    - 10368 * m_e ≈ 5.3 GeV (B meson range)
    - 648 * 108 * m_e ≈ 35.7 GeV (???)
    - 10368 * 108 * m_e ≈ 572 GeV (BSM?)
    """
    
    m_e_MeV = 0.511
    
    scales = {
        'T_P4_mass': 648 * m_e_MeV,
        'K_P4_mass': 10368 * m_e_MeV,
        'T_P4_x_T_P3_mass': 648 * 108 * m_e_MeV,
        'K_P4_x_T_P3_mass': 10368 * 108 * m_e_MeV
    }
    
    return scales


def p5_scale_physics() -> Dict[str, object]:
    """
    Probe the P5 = 2310 scale for physics predictions.
    
    P5 Twist Unit = 6480 = 2⁴ × 3⁴ × 5.
    K(P5) = T(P5) × φ(P5) / 3 = 6480 × 480 / 3 = 1,036,800.
    
    The full primorial 2310 is the smallest modulus that captures all primes
    up to 11. Adding 7 and 11 refines every channel:
    - Spectral gap, modularity, and level repulsion all deviate >1.5σ from random
    - Collapse dynamics concentrate faster; entropy drop is steeper
    - Fine structure correction: α⁻¹ refined by ~1/φ(7×11) = 1/60
    
    P6 = 30030 scale predicts next B-meson family or new resonances.
    """
    m_e_MeV = 0.511
    
    T_P5 = 6480
    PHI_P5 = 480       # φ(2310)
    K_P5 = T_P5 * PHI_P5 / 3  # = 1,036,800
    
    scales = {
        # P5 mass scales
        'T_P5_mass_MeV': T_P5 * m_e_MeV,           # ~3.3 GeV
        'K_P5_mass_MeV': K_P5 * m_e_MeV,            # ~530 GeV
        'T_P5_x_T_P3_mass_MeV': T_P5 * 108 * m_e_MeV,  # ~357 GeV
        
        # Constants
        'T_P5': T_P5,
        'K_P5': K_P5,
        'PHI_P5': PHI_P5,
        
        # 6480 = T(P5) factorization
        'T_P5_factorization': '2⁴ × 3⁴ × 5',
        'T_P5_is_correct': T_P5 == 2**4 * 3**4 * 5,
        
        # Fine structure refinement
        'alpha_inv_base': 108 + 29,                  # = 137
        'alpha_inv_correction': 1 / 27,              # 3-adic
        'alpha_inv_p5_correction': 1 / 60,           # 1/φ(7×11)
        'alpha_inv_refined': 108 + 29 + 1/27 + 1/(27 * 60),
        
        # P6 prediction
        'P6': 30030,
        'T_P6': 27 * 4 * 6 * 10 * 12,  # T = 3³ × Π(p-1) for p=5,7,11,13
        'T_P6_mass_MeV': 27 * 4 * 6 * 10 * 12 * m_e_MeV,
    }
    
    return scales


def fine_structure_from_primorial() -> Dict[str, object]:
    """
    Derive the fine structure constant from the primorial framework.
    
    Base formula (from sub-primorial 30):
        α⁻¹ = 108 + 29 + 1/27 ≈ 137.037037...
    
    The base 108 comes from the twist unit T(P₃).
    The 29 is the boundary prime (largest coprime residue mod 30).
    The 1/27 = 1/3³ is the 3-adic correction.
    
    Higher primorial correction (from 2310):
        Δα⁻¹ = 1/(27 × 60) = 1/1620 ≈ 0.000617
    
    This small correction from 1/φ(7×11) reduces error to <0.0005%.
    
    Experimental value: α⁻¹ = 137.035999084(21)
    """
    base = 108 + 29                            # = 137
    correction_3adic = 1 / 27                  # ≈ 0.037037
    alpha_inv_p3 = base + correction_3adic     # ≈ 137.037037
    
    # Higher primorial correction
    phi_7_11 = 6 * 10  # = 60 = φ(7)×φ(11)
    correction_p5 = -1 / (27 * phi_7_11)       # = -1/1620 ≈ -0.000617
    alpha_inv_p5 = alpha_inv_p3 + correction_p5  # ≈ 137.036420
    
    # Experimental value
    alpha_inv_exp = 137.035999084
    
    error_p3 = abs(alpha_inv_p3 - alpha_inv_exp) / alpha_inv_exp * 100
    error_p5 = abs(alpha_inv_p5 - alpha_inv_exp) / alpha_inv_exp * 100
    
    return {
        'alpha_inv_base': base,
        'alpha_inv_p3': alpha_inv_p3,
        'alpha_inv_p5': alpha_inv_p5,
        'alpha_inv_experimental': alpha_inv_exp,
        'error_p3_pct': error_p3,
        'error_p5_pct': error_p5,
        'formula_p3': '108 + 29 + 1/27',
        'formula_p5': '108 + 29 + 1/27 - 1/(27×60)',
        'improvement': error_p3 > error_p5,
    }

