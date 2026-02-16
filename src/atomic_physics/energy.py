"""
Atomic Energy Calculations.

Implements the Slater Total-Energy Difference method:
    IE = E(ion) - E(neutral) + corrections

Pipeline steps:
    1. Calculate total energy of neutral atom and ion using Slater orbitals
    2. Apply pairing correction for beyond-half-filled subshells (l >= 1)
    3. Apply exchange stabilization for exactly half-filled subshells
    4. Apply relativistic corrections for heavy elements (Z > 36)

Ionization strategy: try removing from ALL subshells, but select
the OUTERMOST positive result. If the Aufbau-last subshell gives
a positive IE, use it. Otherwise fall back to other subshells
ordered by principal quantum number (highest first).
"""

from typing import List, Tuple, Dict, Any
from .constants import (
    RYDBERG_EV, C_PAIR, C_EXCH, C_REL_SP, C_REL_DF,
    N_STAR_VALUES, SHIELDING_COEFFS
)
from .orbitals import electron_configuration
from .shielding import calculate_sigma, get_effective_n


def calculate_total_energy(Z: int, config: List[Tuple[int, int, int]],
                           params: Dict[str, Any] = None) -> float:
    """
    Calculate total electronic energy of a configuration using Slater orbitals.

    E = Σ_i [ n_i × (-R_y × Z_eff_i² / n*_i²) ]
    """
    n_star_vals = params.get('n_star', N_STAR_VALUES) if params else N_STAR_VALUES
    shielding_coeffs = params.get('shielding_coeffs', SHIELDING_COEFFS) if params else SHIELDING_COEFFS

    total_energy = 0.0

    for n, l, count in config:
        sigma = calculate_sigma(n, l, config, coeffs=shielding_coeffs)
        z_eff = max(Z - sigma, 0.01)
        n_star = get_effective_n(n, l, n_star_values=n_star_vals)

        term = -RYDBERG_EV * (z_eff ** 2) / (n_star ** 2)
        total_energy += count * term

    return total_energy


def _compute_ie_for_subshell(Z: int, config_neutral: List[Tuple[int, int, int]],
                              e_neutral: float, idx: int,
                              params: Dict[str, Any] = None) -> float:
    """
    Compute the ionization energy for removing one electron from subshell at idx.
    Returns the IE including all corrections.
    """
    c_pair = params.get('c_pair', C_PAIR) if params else C_PAIR
    c_exch = params.get('c_exch', C_EXCH) if params else C_EXCH
    c_rel_sp = params.get('c_rel_sp', C_REL_SP) if params else C_REL_SP
    c_rel_df = params.get('c_rel_df', C_REL_DF) if params else C_REL_DF

    n, l, count = config_neutral[idx]

    # Create ion configuration
    config_ion = list(config_neutral)
    if count > 1:
        config_ion[idx] = (n, l, count - 1)
    else:
        config_ion = config_ion[:idx] + config_ion[idx + 1:]

    e_ion = calculate_total_energy(Z, config_ion, params)
    ie_base = e_ion - e_neutral

    # ── Step 4: Pairing Correction ──
    # For l >= 1 beyond half-filling, paired electrons experience
    # Coulomb repulsion → easier ionization (lower IE).
    pairing_correction = 0.0
    half_filled = 2 * l + 1
    if l >= 1 and count > half_filled:
        excess = count - half_filled
        pairing_correction = -c_pair * excess / half_filled

    # ── Step 5: Exchange Stabilization ──
    # Exactly half-filled subshells have max exchange energy.
    exchange_correction = 0.0
    if count == half_filled and l >= 1:
        exchange_correction = c_exch * (l / 2.0)

    ie_total = ie_base + pairing_correction + exchange_correction

    # ── Step 6: Relativistic Correction (Z > 36) ──
    if Z > 36:
        alpha = 1.0 / 137.036
        rel_term = Z ** 2 * alpha ** 2
        if l <= 1:
            ie_total *= (1.0 + c_rel_sp * rel_term)
        else:
            ie_total *= (1.0 - c_rel_df * rel_term)

    return ie_total


def calculate_ionization_transition(Z: int, config_start: List[Tuple[int, int, int]],
                                  params: Dict[str, Any] = None) -> Tuple[float, List[Tuple[int, int, int]]]:
    """
    Calculate the ionization energy and resulting configuration for the next ionization step.
    
    Returns:
        (ie, config_ion)
    """
    if not config_start:
        return 0.0, []

    e_start = calculate_total_energy(Z, config_start, params)

    # Compute IE for all subshells
    candidates = []  # (ie, n, l, idx, config_ion)
    
    for idx, (n, l, count) in enumerate(config_start):
        # Create ion configuration
        config_ion = list(config_start)
        if count > 1:
            config_ion[idx] = (n, l, count - 1)
        else:
            config_ion = config_ion[:idx] + config_ion[idx + 1:]
            
        ie = _compute_ie_for_subshell(Z, config_start, e_start, idx, params)
        
        if ie > 0:
            candidates.append((ie, n, l, idx, config_ion))

    if not candidates:
        # Should not happen for neutral atoms, but possible for highly stripped ions
        return 0.0, []

    # Prefer outermost subshell: sort by highest n first, then lowest IE
    last_idx = len(config_start) - 1
    n_last, l_last, _ = config_start[last_idx]

    # Try valence candidates: last subshell and the one before if d/f
    valence_candidates = []
    for ie, n, l, idx, cfg in candidates:
        if idx == last_idx:
            valence_candidates.append((ie, n, l, cfg))
        elif idx == last_idx - 1 and l_last >= 2:
            valence_candidates.append((ie, n, l, cfg))

    if valence_candidates:
        # If we have both an s orbital and a d/f orbital as candidates,
        # prefer the s orbital (it's the one that ionizes first in TMs)
        s_candidates = [(ie, cfg) for ie, n, l, cfg in valence_candidates if l == 0]
        if s_candidates and l_last >= 2:
             # Find min IE among s candidates
             return min(s_candidates, key=lambda x: x[0])
        
        # Otherwise min IE among all valence candidates
        best = min(valence_candidates, key=lambda x: x[0])
        return best[0], best[3]

    # Fallback: use minimum positive IE from any subshell
    best = min(candidates, key=lambda x: x[0])
    return best[0], best[4]


def ionization_energy(Z: int, params: Dict[str, Any] = None) -> float:
    """
    Predict first ionization energy using the Slater Total-Energy Difference method.
    """
    config_neutral = electron_configuration(Z)
    if not config_neutral:
        return 0.0
        
    ie, _ = calculate_ionization_transition(Z, config_neutral, params)
    return ie


def successive_ionization_energies(Z: int, max_k: int = 10, 
                                 params: Dict[str, Any] = None) -> List[float]:
    """
    Calculate successive ionization energies (IE_1, IE_2, ..., IE_k).
    
    Args:
        Z: Atomic number
        max_k: Maximum ionization level to compute (default 10)
        params: Optional model parameters
        
    Returns:
        List of ionization energies [IE_1, IE_2, ...]
    """
    ies = []
    current_config = electron_configuration(Z)
    
    # Determine how many electrons we can remove
    total_electrons = sum(c for _, _, c in current_config)
    steps = min(max_k, total_electrons)
    
    for _ in range(steps):
        ie, next_config = calculate_ionization_transition(Z, current_config, params)
        if ie <= 0:
            break
        ies.append(ie)
        current_config = next_config
        
    return ies
