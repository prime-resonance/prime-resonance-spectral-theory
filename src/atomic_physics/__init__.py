"""
Atomic Physics Package.

Implements the Justified Slater ΔE Pipeline for predicting
ionization energies from first principles with physically
motivated corrections.

Modules:
    constants   - Tunable global parameters (n*, shielding, corrections)
    orbitals    - Electron configuration and orbital structure
    shielding   - Slater shielding / screening calculations
    energy      - Total energy and ionization energy calculations
    optimization - Grid search parameter optimization
    reporting   - Periodic table generation and accuracy statistics
"""

from .constants import RYDBERG_EV
from .orbitals import (
    orbital_quantum_numbers,
    electron_configuration,
    principal_quantum_number,
    validate_orbital_structure,
)
from .shielding import (
    get_slater_group_index,
    get_effective_n,
    calculate_sigma,
)
from .energy import (
    calculate_total_energy,
    ionization_energy,
    calculate_ionization_transition,
    successive_ionization_energies,
)
from .reporting import (
    generate_periodic_table,
    generate_periodic_table_summary,
    alkali_ionization_trend,
    config_to_string,
)

__all__ = [
    'RYDBERG_EV',
    'orbital_quantum_numbers',
    'electron_configuration',
    'principal_quantum_number',
    'validate_orbital_structure',
    'get_slater_group_index',
    'get_effective_n',
    'calculate_sigma',
    'calculate_total_energy',
    'ionization_energy',
    'calculate_ionization_transition',
    'successive_ionization_energies',
    'generate_periodic_table',
    'generate_periodic_table_summary',
    'alkali_ionization_trend',
    'config_to_string',
]
