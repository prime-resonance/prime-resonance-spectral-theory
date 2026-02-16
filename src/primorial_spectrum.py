"""
Primorial-Aligned Spectral Analysis.

Tests the key prediction: when the Legendre Hamiltonian's basis size matches
a primorial-aligned count (φ(90)=24, φ(630)=144, etc.), the trace invariants
and spectral properties should encode the corresponding twist unit.

Also explores the relationship between the twist unit hierarchy
{27, 108, 648, 6480} and the 288-element symmetry of the figure-eight knot.
"""

from typing import Dict, List, Tuple
import numpy as np
from .prime_utils import sieve_primes, legendre_symbol, digital_root
from .collapse_hamiltonian import build_full_hamiltonian
from .analytic_108 import euler_totient, general_twist_formula


# Primorial-aligned basis sizes
PRIMORIAL_ALIGNMENTS = {
    'P2': {'primorial': 6, 'period_mod': 18, 'phi': 6, 'twist': 27},
    'P3': {'primorial': 30, 'period_mod': 90, 'phi': 24, 'twist': 108},
    'P4': {'primorial': 210, 'period_mod': 630, 'phi': 144, 'twist': 648},
    'P5': {'primorial': 2310, 'period_mod': 6930, 'phi': 1440, 'twist': 6480},
}


def primorial_aligned_traces(
    max_prime: int = 1000,
    gamma: float = 1.0,
) -> Dict[str, Dict[str, object]]:
    """
    Compute Hamiltonian trace invariants at primorial-aligned basis sizes.

    For each primorial P_k, use the first φ(lcm(P_k, 9)) odd primes as the basis.
    Compute Tr(H), Tr(H²), Tr(H³) and look for the twist unit T(P_k).

    Returns:
        Dict of {primorial_name: {basis_size, traces, twist_unit, ...}}
    """
    all_primes = [p for p in sieve_primes(max_prime) if p > 2]

    results = {}
    for name, info in PRIMORIAL_ALIGNMENTS.items():
        n = info['phi']
        if n > len(all_primes):
            continue

        primes = all_primes[:n]
        H = build_full_hamiltonian(primes, gamma, use_legendre=True)

        # Compute traces
        traces = {}
        H_power = np.eye(n, dtype=np.complex128)
        for k in range(1, 5):
            H_power = H_power @ H
            traces[k] = complex(np.trace(H_power))

        # Eigenvalues
        eigs = np.linalg.eigvals(H)
        real_eigs = np.sort(eigs.real)

        # Spectral width
        width = float(real_eigs[-1] - real_eigs[0])

        # Digital root of spectral width (mod twist unit)
        width_mod_twist = width % info['twist'] if info['twist'] > 0 else 0

        results[name] = {
            'basis_size': n,
            'primes': primes,
            'twist_unit': info['twist'],
            'traces': traces,
            'spectral_width': width,
            'width_mod_twist': width_mod_twist,
            'eigenvalues_real': real_eigs.tolist(),
            'n_eigenvalues': len(eigs),
            'trace_1_abs': abs(traces[1]),
            'trace_2_abs': abs(traces[2]),
        }

    return results


def twist_unit_in_trace_ratios(
    max_prime: int = 1000,
    gamma: float = 1.0,
) -> Dict[str, object]:
    """
    Search for twist units in the RATIOS of trace invariants
    across different primorial-aligned basis sizes.

    The prediction: Tr_P4(H²) / Tr_P3(H²) ≈ T(P4)/T(P3) = 648/108 = 6.

    Returns:
        Dict with trace ratios and their proximity to twist unit ratios.
    """
    traces = primorial_aligned_traces(max_prime, gamma)

    ratios = {}
    names = sorted(traces.keys())

    for i in range(len(names) - 1):
        n1, n2 = names[i], names[i + 1]
        t1, t2 = traces[n1], traces[n2]

        for k in [1, 2, 3]:
            if abs(t1['traces'][k]) > 1e-10:
                ratio = abs(t2['traces'][k]) / abs(t1['traces'][k])
                twist_ratio = t2['twist_unit'] / t1['twist_unit']
                ratios[f'Tr{k}_{n2}/{n1}'] = {
                    'trace_ratio': ratio,
                    'twist_ratio': twist_ratio,
                    'relative_error': abs(ratio - twist_ratio) / twist_ratio if twist_ratio != 0 else float('inf'),
                }

    return ratios


def the_288_connection() -> Dict[str, object]:
    """
    Explore the relationship between T(P₃)=108 and the 288-element symmetry.

    From the formalism: 288 = s×c×b×u × 8 = 36 × 8.
    Trefoil invariants: s=6, c=3, b=2, u=1 → product = 36.

    Key relationships:
    - 288/108 = 8/3 (from formalism)
    - 648/288 = 9/4 = (3/2)²
    - 108 × φ(30)/3 = 108 × 8/3 = 288

    Tests if the 288-element symmetry has a clean primorial-level interpretation.
    """
    # Basic constants
    T_P3 = 108  # = 2² × 3³
    sym_288 = 288  # = 2⁵ × 3² (figure-eight knot symmetry)
    T_P4 = 648  # = 2³ × 3⁴

    # The trefoil invariant product
    s, c, b, u = 6, 3, 2, 1
    trefoil_product = s * c * b * u  # = 36
    assert trefoil_product == 36
    assert sym_288 == trefoil_product * 8

    # Key ratios
    ratio_288_108 = sym_288 / T_P3  # = 8/3
    ratio_648_288 = T_P4 / sym_288  # = 9/4
    ratio_648_108 = T_P4 / T_P3    # = 6 = φ(7)

    # Relationship: 288 = T(P₃) × φ(P₃) / 3
    phi_30 = euler_totient(30)  # = 8
    assert T_P3 * phi_30 / 3 == sym_288  # 108 × 8/3 = 288

    # Generalization: K(P_k) = T(P_k) × φ(P_k) / 3
    K_P3 = T_P3 * phi_30 / 3  # = 288
    phi_210 = euler_totient(210)  # = 48
    K_P4 = T_P4 * phi_210 / 3  # = 648 × 48 / 3 = 10368

    return {
        'T_P3': T_P3,
        'sym_288': sym_288,
        'T_P4': T_P4,
        'trefoil_product': trefoil_product,
        'ratio_288_108': ratio_288_108,
        'ratio_648_288': ratio_648_288,
        'ratio_648_108': ratio_648_108,
        'phi_30': phi_30,
        'relationship': f'288 = {T_P3} × {phi_30} / 3',
        'K_P3': K_P3,
        'K_P4': K_P4,
        'K_P4_factorization': f'{int(K_P4)} = 2^7 × 3^4' if K_P4 == 2**7 * 3**4 else f'{int(K_P4)}',
    }


def mass_ratio_search(twist_unit: int, max_multiplier: int = 50) -> List[Dict[str, object]]:
    """
    Search for known particle mass ratios (to electron) that are approximately
    n × twist_unit ± correction, where the correction involves small primes.

    This tests whether the P₄ twist unit (648) predicts known particle masses
    that are currently "unexplained" by the P₃ formulas.

    Returns:
        List of candidate mass matches.
    """
    # Known mass ratios to electron
    particles = {
        'electron': 1.0,
        'muon': 206.768,
        'tau': 3477.23,
        'proton': 1836.153,
        'neutron': 1838.684,
        'W_boson': 157296.0,  # 80.4 GeV / 0.511 MeV
        'Z_boson': 178449.0,  # 91.2 GeV / 0.511 MeV
        'higgs': 244912.0,    # 125.1 GeV / 0.511 MeV
        'charm_quark': 2494.0,  # ~1275 MeV
        'bottom_quark': 8180.0,  # ~4180 MeV
        'top_quark': 338646.0,  # 173 GeV
        'pion_charged': 273.1,  # 139.6 MeV
        'pion_neutral': 264.1,  # 135.0 MeV
        'kaon': 966.0,  # 493.7 MeV
    }

    candidates = []
    for name, ratio in particles.items():
        n = round(ratio / twist_unit)
        if n == 0:
            continue
        predicted = n * twist_unit
        residual = ratio - predicted
        error_pct = abs(residual) / ratio * 100

        # Also check n × twist_unit ± small prime corrections
        best_correction = None
        best_error = error_pct
        correction_type = "none"

        # Standard cubic/quadratic corrections
        standard_corrections = [0, 1, -1, 3, -3, 4, -4, 9, -9, 27, -27]
        for correction in standard_corrections:
            candidate = n * twist_unit + correction
            err = abs(ratio - candidate) / ratio * 100
            if err < best_error:
                best_error = err
                best_correction = correction
                correction_type = "standard"

        # Special composite correction for heavy leptons (Tau)
        # Combines Quintic (5^3) and Quadratic (2^2) channels
        # Formula: n * 108 + 125 + 4
        quintic_quad = 125 + 4
        candidate_q = (n - 1) * twist_unit + quintic_quad # Adjust n because 125+4 > 108
        err_q = abs(ratio - candidate_q) / ratio * 100
        
        # Also check without n-1 adjustment (if residue is small)
        candidate_q2 = n * twist_unit + quintic_quad
        err_q2 = abs(ratio - candidate_q2) / ratio * 100

        if err_q < best_error:
            best_error = err_q
            best_correction = quintic_quad
            correction_type = "quintic_quadratic_shifted" # n -> n-1
        
        if err_q2 < best_error:
            best_error = err_q2
            best_correction = quintic_quad
            correction_type = "quintic_quadratic"

        candidates.append({
            'particle': name,
            'mass_ratio': ratio,
            'n': n,
            'n_times_twist': n * twist_unit,
            'residual': residual,
            'error_pct': error_pct,
            'best_correction': best_correction,
            'correction_type': correction_type,
            'best_error_pct': best_error,
        })

    return sorted(candidates, key=lambda c: c['best_error_pct'])
