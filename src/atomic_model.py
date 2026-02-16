"""
Atomic Model Generator from the Primorial Reciprocity Framework.

Generates the periodic table of elements — orbital structure, quantum numbers,
ionization energies — using the Justified Slater ΔE Pipeline.

Core principle: The same primorial decomposition 2310 = 2 × 3 × 5 × 7 × 11
that governs particle physics also governs atomic structure.

Pipeline:
    1. Base energy via Slater Total-Energy Difference: IE = E(ion) - E(neutral)
    2. Effective quantum numbers n*(n,l) with orbital penetration effects
    3. l-dependent shielding constants (s penetrates core more than p)
    4. Pairing correction for > half-filled subshells
    5. Exchange stabilization for exactly half-filled subshells
    6. Relativistic corrections for Z > 36

This module is a facade that delegates to the atomic_physics subpackage.
"""

from typing import Dict, List, Tuple, Optional

# Re-export constants from the original model (kept for backward compat)
from .atomic_physics.constants import RYDBERG_EV

# Primorial framework constants (informational)
T_P3 = 108
T_P5 = 6480
PHI_P3 = 8       # φ(30) = 8 coprime residue classes
PHI_P5 = 480     # φ(2310) = 480

# Screening denominators from the primorial sieve
SCREENING_DENOMS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15,
                    16, 18, 20, 21, 22, 24, 25, 27, 28, 30, 33, 35,
                    36, 40, 42, 44, 45, 48, 49, 50, 54, 55, 56, 60,
                    63, 66, 70, 72, 75, 77, 80, 84, 88, 90, 98, 99,
                    100, 105, 108, 110, 120, 125, 126, 132, 140, 144,
                    150, 154, 160, 165, 168, 175, 176, 180, 189, 196,
                    198, 200, 210, 220, 225, 231, 240, 250, 252, 264,
                    270, 275, 280, 288, 294, 300, 308, 315, 330, 336,
                    350, 360, 375, 378, 385, 396, 400, 420, 432, 440,
                    450, 462, 480, 490, 495, 500, 504, 525, 528, 540,
                    550, 560, 567, 576, 588, 600, 616, 630, 648, 660,
                    675, 693, 700, 720, 750, 756, 770, 792, 800, 840,
                    864, 880, 900, 924, 945, 960, 980, 990, 1000, 1008,
                    1050, 1056, 1080, 1100, 1120, 1125, 1134, 1155, 1176,
                    1200, 1232, 1260, 1296, 1320, 1350, 1386, 1400, 1440,
                    1500, 1512, 1540, 1575, 1584, 1600, 1620, 1680, 1728,
                    1760, 1764, 1800, 1848, 1890, 1920, 1925, 1960, 1980,
                    2000, 2016, 2100, 2160, 2200, 2250, 2268, 2310]

# Import NIST data
from .nist_data import (
    NIST_IONIZATION_ENERGIES,
    SHELL_FILLING_ORDER,
    NOBLE_GAS_Z,
    ALKALI_Z,
    PERIOD_LENGTHS,
    get_ionization_energy,
    get_element_symbol,
    get_period,
)

# Import all pipeline functions from the subpackage
from .atomic_physics.orbitals import (
    orbital_quantum_numbers,
    electron_configuration,
    principal_quantum_number,
    validate_orbital_structure,
)
from .atomic_physics.shielding import (
    get_slater_group_index as _slater_group_index,
    get_effective_n as _get_effective_n,
    get_shielding_coefficient as _get_shielding_coefficient,
    calculate_sigma,
)
from .atomic_physics.energy import (
    calculate_total_energy,
    ionization_energy as _ionization_energy,
    calculate_ionization_transition,
    successive_ionization_energies,
)
from .atomic_physics.reporting import (
    generate_periodic_table,
    generate_periodic_table_summary,
    alkali_ionization_trend,
    config_to_string,
)


def z_effective_primorial(Z: int) -> float:
    """
    Compute Z_eff for the outermost electron (for display purposes).
    """
    config = electron_configuration(Z)
    if not config:
        return 0.0
    n, l, count = config[-1]
    sigma = calculate_sigma(n, l, config)
    return Z - sigma


def ionization_energy_primorial(Z: int) -> float:
    """
    Predict ionization energy using the Justified Slater ΔE Pipeline.
    
    Delegates to atomic_physics.energy.ionization_energy().
    """
    return _ionization_energy(Z)
