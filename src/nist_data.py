"""
NIST Atomic Ionization Energy Data.

Reference ionization energies (ground state, first ionization) for all
elements Z=1..118 from the NIST Atomic Spectra Database.

Source: NIST ASD (https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html)
Units: eV
Last updated: 2024

These values are used to validate the primorial reciprocity framework's
parameter-free predictions for atomic physics.
"""

from typing import Dict, Tuple, Optional

# NIST first ionization energies in eV for Z=1..118
# Format: Z -> (symbol, name, ionization_energy_eV)
NIST_IONIZATION_ENERGIES: Dict[int, Tuple[str, str, float]] = {
    1:  ('H',  'Hydrogen',      13.598434),
    2:  ('He', 'Helium',        24.587388),
    3:  ('Li', 'Lithium',        5.391715),
    4:  ('Be', 'Beryllium',      9.322700),
    5:  ('B',  'Boron',          8.298020),
    6:  ('C',  'Carbon',        11.260300),
    7:  ('N',  'Nitrogen',      14.534140),
    8:  ('O',  'Oxygen',        13.618055),
    9:  ('F',  'Fluorine',      17.422820),
    10: ('Ne', 'Neon',          21.564540),
    11: ('Na', 'Sodium',         5.139076),
    12: ('Mg', 'Magnesium',      7.646236),
    13: ('Al', 'Aluminium',      5.985768),
    14: ('Si', 'Silicon',        8.151683),
    15: ('P',  'Phosphorus',    10.486686),
    16: ('S',  'Sulfur',        10.360010),
    17: ('Cl', 'Chlorine',      12.967632),
    18: ('Ar', 'Argon',         15.759610),
    19: ('K',  'Potassium',      4.340664),
    20: ('Ca', 'Calcium',        6.113155),
    21: ('Sc', 'Scandium',       6.561490),
    22: ('Ti', 'Titanium',       6.828120),
    23: ('V',  'Vanadium',       6.746187),
    24: ('Cr', 'Chromium',       6.766510),
    25: ('Mn', 'Manganese',      7.434018),
    26: ('Fe', 'Iron',           7.902468),
    27: ('Co', 'Cobalt',         7.881010),
    28: ('Ni', 'Nickel',         7.639877),
    29: ('Cu', 'Copper',         7.726380),
    30: ('Zn', 'Zinc',           9.394199),
    31: ('Ga', 'Gallium',        5.999302),
    32: ('Ge', 'Germanium',      7.899435),
    33: ('As', 'Arsenic',        9.788600),
    34: ('Se', 'Selenium',       9.752380),
    35: ('Br', 'Bromine',       11.813800),
    36: ('Kr', 'Krypton',       13.999610),
    37: ('Rb', 'Rubidium',       4.177128),
    38: ('Sr', 'Strontium',      5.694867),
    39: ('Y',  'Yttrium',        6.217300),
    40: ('Zr', 'Zirconium',      6.634050),
    41: ('Nb', 'Niobium',        6.758850),
    42: ('Mo', 'Molybdenum',     7.092430),
    43: ('Tc', 'Technetium',     7.119380),
    44: ('Ru', 'Ruthenium',      7.360500),
    45: ('Rh', 'Rhodium',        7.458900),
    46: ('Pd', 'Palladium',      8.336900),
    47: ('Ag', 'Silver',         7.576234),
    48: ('Cd', 'Cadmium',        8.993820),
    49: ('In', 'Indium',         5.786360),
    50: ('Sn', 'Tin',            7.343918),
    51: ('Sb', 'Antimony',       8.608400),
    52: ('Te', 'Tellurium',      9.009660),
    53: ('I',  'Iodine',        10.451260),
    54: ('Xe', 'Xenon',         12.129843),
    55: ('Cs', 'Cesium',         3.893905),
    56: ('Ba', 'Barium',         5.211664),
    57: ('La', 'Lanthanum',      5.576900),
    58: ('Ce', 'Cerium',         5.538600),
    59: ('Pr', 'Praseodymium',   5.473000),
    60: ('Nd', 'Neodymium',      5.525000),
    61: ('Pm', 'Promethium',     5.582000),
    62: ('Sm', 'Samarium',       5.643700),
    63: ('Eu', 'Europium',       5.670385),
    64: ('Gd', 'Gadolinium',     6.149800),
    65: ('Tb', 'Terbium',        5.863700),
    66: ('Dy', 'Dysprosium',     5.938900),
    67: ('Ho', 'Holmium',        6.021500),
    68: ('Er', 'Erbium',         6.107300),
    69: ('Tm', 'Thulium',        6.184310),
    70: ('Yb', 'Ytterbium',      6.254160),
    71: ('Lu', 'Lutetium',       5.425871),
    72: ('Hf', 'Hafnium',        6.825070),
    73: ('Ta', 'Tantalum',       7.549570),
    74: ('W',  'Tungsten',       7.864050),
    75: ('Re', 'Rhenium',        7.833520),
    76: ('Os', 'Osmium',         8.438200),
    77: ('Ir', 'Iridium',        8.967000),
    78: ('Pt', 'Platinum',       8.958830),
    79: ('Au', 'Gold',           9.225554),
    80: ('Hg', 'Mercury',       10.437505),
    81: ('Tl', 'Thallium',       6.108194),
    82: ('Pb', 'Lead',           7.416600),
    83: ('Bi', 'Bismuth',        7.285500),
    84: ('Po', 'Polonium',       8.414000),
    85: ('At', 'Astatine',       9.317500),
    86: ('Rn', 'Radon',         10.748500),
    87: ('Fr', 'Francium',       4.072742),
    88: ('Ra', 'Radium',         5.278423),
    89: ('Ac', 'Actinium',       5.380226),
    90: ('Th', 'Thorium',        6.306600),
    91: ('Pa', 'Protactinium',   5.890000),
    92: ('U',  'Uranium',        6.194050),
    93: ('Np', 'Neptunium',      6.265700),
    94: ('Pu', 'Plutonium',      6.026000),
    95: ('Am', 'Americium',      5.973800),
    96: ('Cm', 'Curium',         5.991400),
    97: ('Bk', 'Berkelium',      6.197900),
    98: ('Cf', 'Californium',    6.281600),
    99: ('Es', 'Einsteinium',    6.367400),
    100: ('Fm', 'Fermium',       6.500000),
    101: ('Md', 'Mendelevium',   6.580000),
    102: ('No', 'Nobelium',      6.626200),
    103: ('Lr', 'Lawrencium',    4.960100),
    104: ('Rf', 'Rutherfordium', 6.010000),
    105: ('Db', 'Dubnium',       6.890000),
    106: ('Sg', 'Seaborgium',    7.100000),
    107: ('Bh', 'Bohrium',       7.700000),
    108: ('Hs', 'Hassium',       7.600000),
    109: ('Mt', 'Meitnerium',    9.100000),
    110: ('Ds', 'Darmstadtium',  9.900000),
    111: ('Rg', 'Roentgenium',  10.600000),
    112: ('Cn', 'Copernicium',  11.970000),
    113: ('Nh', 'Nihonium',      7.306000),
    114: ('Fl', 'Flerovium',     8.539000),
    115: ('Mc', 'Moscovium',     5.580000),
    116: ('Lv', 'Livermorium',   6.780000),
    117: ('Ts', 'Tennessine',    7.700000),
    118: ('Og', 'Oganesson',     8.914000),
}

# Electron configuration shells: n, l, max_electrons
# Standard quantum numbers for ground-state filling order
SHELL_FILLING_ORDER = [
    (1, 0, 2),   # 1s
    (2, 0, 2),   # 2s
    (2, 1, 6),   # 2p
    (3, 0, 2),   # 3s
    (3, 1, 6),   # 3p
    (4, 0, 2),   # 4s
    (3, 2, 10),  # 3d
    (4, 1, 6),   # 4p
    (5, 0, 2),   # 5s
    (4, 2, 10),  # 4d
    (5, 1, 6),   # 5p
    (6, 0, 2),   # 6s
    (4, 3, 14),  # 4f
    (5, 2, 10),  # 5d
    (6, 1, 6),   # 6p
    (7, 0, 2),   # 7s
    (5, 3, 14),  # 5f
    (6, 2, 10),  # 6d
    (7, 1, 6),   # 7p
]

# Noble gas Z values (period closings)
NOBLE_GAS_Z = [2, 10, 18, 36, 54, 86, 118]

# Alkali metal Z values (period openings, single valence electron)
ALKALI_Z = [1, 3, 11, 19, 37, 55, 87]

# Period lengths (number of elements per period)
PERIOD_LENGTHS = [2, 8, 8, 18, 18, 32, 32]


def get_ionization_energy(Z: int) -> Optional[float]:
    """Get NIST ionization energy for element Z."""
    if Z in NIST_IONIZATION_ENERGIES:
        return NIST_IONIZATION_ENERGIES[Z][2]
    return None


def get_element_symbol(Z: int) -> Optional[str]:
    """Get element symbol for atomic number Z."""
    if Z in NIST_IONIZATION_ENERGIES:
        return NIST_IONIZATION_ENERGIES[Z][0]
    return None


def get_element_name(Z: int) -> Optional[str]:
    """Get element name for atomic number Z."""
    if Z in NIST_IONIZATION_ENERGIES:
        return NIST_IONIZATION_ENERGIES[Z][1]
    return None


def get_period(Z: int) -> int:
    """Get the period number for element Z."""
    cumulative = 0
    for period, length in enumerate(PERIOD_LENGTHS, 1):
        cumulative += length
        if Z <= cumulative:
            return period
    return 7  # Default to period 7 for superheavy


def get_group_type(Z: int) -> str:
    """Classify element by type based on Z."""
    if Z in NOBLE_GAS_Z:
        return 'noble_gas'
    if Z in ALKALI_Z:
        return 'alkali'
    if Z in [4, 12, 20, 38, 56, 88]:
        return 'alkaline_earth'
    if Z in range(57, 72) or Z in range(89, 104):
        return 'lanthanide_actinide'
    if Z in range(21, 31) or Z in range(39, 49) or Z in range(72, 81):
        return 'transition_metal'
    return 'main_group'
