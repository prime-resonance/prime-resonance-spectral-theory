"""
Slater Shielding and Effective Quantum Numbers.

Implements l-dependent shielding with source-orbital awareness:
- s-orbitals penetrate the core, so inner electrons shield them LESS (0.85)
- d/f source orbitals are compact, so they shield outer s LESS than s/p sources do
- Standard Slater grouping: (1s)(2s,2p)(3s,3p)(3d)(4s,4p)(4d)(4f)...
- 3-adic tower correction: valence s-orbitals (n > 2) receive a shielding
  reduction of -1/3^(n-2), reflecting the cubic reciprocity channel's
  role in penetration depth across successive shells.

Key physical insight: The shielding effectiveness depends on BOTH the
target orbital's penetration AND the source orbital's spatial extent.
"""

from typing import List, Tuple, Dict, Any
from .constants import N_STAR_VALUES, SHIELDING_COEFFS


def get_slater_group_index(n: int, l: int) -> int:
    """
    Assign a Slater group index for orbital (n, l).

    Slater groups: (1s)(2s,2p)(3s,3p)(3d)(4s,4p)(4d)(4f)...
    s,p orbitals with same n share a group (index = n*10).
    d,f orbitals each get their own group (index = n*10 + l).
    
    This ordering ensures 3d (32) < 4s (40) < 4d (42) < 4f (43) < 5s (50).
    """
    if l <= 1:
        return n * 10  # e.g. 1s->10, 2s/2p->20, 3s/3p->30
    else:
        return n * 10 + l  # e.g. 3d->32, 4f->43


def get_effective_n(n: int, l: int,
                    n_star_values: Dict[int, Dict[int, float]] = None) -> float:
    """
    Get effective quantum number n* based on (n, l).

    Physical justification for n*(n,l):
    - s-orbitals penetrate core → lower n* → tighter binding
    - d,f orbitals are compact → lower n* than n
    - 6s has relativistic contraction → reduced n*
    """
    if n_star_values is None:
        n_star_values = N_STAR_VALUES

    if n in n_star_values and l in n_star_values[n]:
        return n_star_values[n][l]

    # Fallback for quantum numbers not in the table
    return float(n)


def get_shielding_coefficient(target_l: int, source_l: int, relation: str,
                              coeffs: Dict[str, Dict[str, float]] = None) -> float:
    """
    Get shielding coefficient based on target l, source l, and relationship.

    The coefficient c(g_j, g_i, l_i) depends on:
    - target_l: angular momentum of the electron being shielded
    - source_l: angular momentum of the shielding electron
    - relation: 'same_group', 'one_below', 'deep_core'

    Physical justification:
    - s penetrates core → inner electrons shield it LESS (0.85 vs 0.93 for p)
    - d/f sources are compact → they shield outer s/p LESS than s/p sources
    """
    if coeffs is None:
        coeffs = SHIELDING_COEFFS

    # Determine target type
    if target_l == 0:
        l_type = 's'
    elif target_l == 1:
        l_type = 'p'
    else:
        l_type = 'df'

    base_coeff = coeffs[l_type].get(relation, 1.00)

    # Source-l refinement: d/f sources shield s/p targets less effectively
    # in the "one_below" relationship because d/f orbitals are compact
    # and don't overlap as much with the penetrating s/p orbitals.
    if relation == 'one_below' and source_l >= 2 and target_l <= 1:
        # Apply a reduction factor for d/f → s/p shielding
        # This captures the fact that 3d doesn't shield 4s as well as 3s/3p do
        df_reduction = coeffs.get('df_source_reduction', {}).get(l_type, 1.0)
        base_coeff *= df_reduction

    return base_coeff


def calculate_sigma(target_n: int, target_l: int,
                    config: List[Tuple[int, int, int]],
                    coeffs: Dict[str, Dict[str, float]] = None) -> float:
    """
    Calculate Slater shielding constant sigma for a specific subshell.

    Uses Slater grouping with l-dependent coefficients:
    1. Electrons in higher groups: 0 contribution
    2. Electrons in same group: same_group coefficient
    3. For s/p targets: one_below = (n-1) shell, deep_core = (n-2) and below
    4. For d/f targets: all lower groups contribute ~1.00

    The (n-1) shell for s/p targets means ALL electrons with
    principal quantum number n-1 (regardless of their l), but
    d/f sources get reduced effectiveness via get_shielding_coefficient.
    """
    if coeffs is None:
        coeffs = SHIELDING_COEFFS

    target_group = get_slater_group_index(target_n, target_l)
    sigma = 0.0

    for n, l, count in config:
        group = get_slater_group_index(n, l)

        # Self-interaction correction: exclude one electron from target subshell
        if n == target_n and l == target_l:
            count_contribution = count - 1
        else:
            count_contribution = count

        if count_contribution <= 0:
            continue

        # Electrons in higher groups contribute 0
        if group > target_group:
            continue

        # Determine the relationship
        if group == target_group:
            relation = 'same_group'
        elif target_l <= 1:
            # For s/p targets: (n-1) shell = one_below, (n-2) and lower = deep_core
            # "Shell" means principal quantum number here
            if n == target_n - 1:
                relation = 'one_below'
            else:
                relation = 'deep_core'
        else:
            # For d/f targets: all lower groups contribute ~1.00
            # We still distinguish for flexibility, but default coefficients are 1.00
            if n >= target_n - 1:
                relation = 'one_below'
            else:
                relation = 'deep_core'

        coeff = get_shielding_coefficient(target_l, l, relation, coeffs)
        sigma += count_contribution * coeff

    # ── 3-adic tower correction for valence s-orbitals ──
    # Physical basis: The cubic reciprocity channel (prime 3) governs
    # how deeply an s-orbital penetrates successive shells. For n > 2,
    # each additional shell adds a correction of -1/3^(n-2) to σ,
    # reflecting the 3-adic hierarchy of shielding depth.
    # This is a structurally forced correction (zero new free parameters).
    if target_l == 0 and target_n > 2:
        tower_height = target_n - 2
        delta_sigma = -1.0 / (3 ** tower_height)
        sigma += delta_sigma

    # ── 2-adic valence-s penetration term ──
    # Analogous to the 3-adic but for prime-2 channel.
    # Small fixed prefactor 0.075, tower height n-2.
    if target_l == 0 and target_n >= 3:
        delta_sigma_2 = -0.075 / (2 ** (target_n - 2))
        sigma += delta_sigma_2

    return sigma
