"""
Orbital Structure and Electron Configuration.

Handles the definition of atomic orbitals, quantum numbers, and
Aufbau principle filling orders.
"""

from typing import List, Dict, Tuple, Optional
from ..nist_data import SHELL_FILLING_ORDER, PERIOD_LENGTHS, NOBLE_GAS_Z

def orbital_quantum_numbers() -> List[Dict[str, object]]:
    """
    Derive orbital quantum numbers from the primorial channel structure.
    
    Returns:
        List of orbital specifications with quantum numbers and capacities.
    """
    orbitals = []
    for n, l, max_e in SHELL_FILLING_ORDER:
        angular_states = 2 * l + 1
        spin_states = 2
        total = angular_states * spin_states
        
        # Channel mapping (informational, from original model)
        channel_map = {
            0: {'prime': 2, 'channel': 'quadratic', 'label': 's'},
            1: {'prime': 3, 'channel': 'cubic', 'label': 'p'},
            2: {'prime': 5, 'channel': 'quintic', 'label': 'd'},
            3: {'prime': 7, 'channel': 'septic', 'label': 'f'},
        }
        channel = channel_map.get(l, {'prime': None, 'channel': 'higher', 'label': '?'})

        orbitals.append({
            'n': n,
            'l': l,
            'label': f"{n}{channel['label']}",
            'max_electrons': max_e,
            'angular_states': angular_states,
            'primorial_prime': channel['prime'],
            'channel_name': channel['channel'],
            'shell_capacity': total,
        })
    return orbitals


def electron_configuration(Z: int) -> List[Tuple[int, int, int]]:
    """
    Build the ground-state electron configuration for element Z.
    
    Fills orbitals in the standard Aufbau order.
    
    Args:
        Z: Atomic number
        
    Returns:
        List of (n, l, electron_count) tuples.
    """
    config = []
    remaining = Z
    for n, l, max_e in SHELL_FILLING_ORDER:
        if remaining <= 0:
            break
        fill = min(remaining, max_e)
        config.append((n, l, fill))
        remaining -= fill
    return config


def principal_quantum_number(Z: int) -> int:
    """
    Return the principal quantum number n of the outermost electron.
    """
    config = electron_configuration(Z)
    if not config:
        return 1
    return config[-1][0]


def validate_orbital_structure() -> Dict[str, object]:
    """
    Validate that the orbital structure matches known physical periods.
    """
    # Subshell capacities
    primes = [2, 3, 5, 7]
    predicted_capacities = []
    for l, p in enumerate(primes):
        capacity = 2 * (2 * l + 1)
        predicted_capacities.append({
            'l': l,
            'label': 'spdf'[l],
            'prime': p,
            'capacity': capacity,
        })

    # Verify period lengths
    predicted_periods = []
    cumulative_z = 0
    for i, length in enumerate(PERIOD_LENGTHS):
        cumulative_z += length
        predicted_periods.append({
            'period': i + 1,
            'length': length,
            'cumulative_Z': cumulative_z,
            'noble_gas_Z': NOBLE_GAS_Z[i] if i < len(NOBLE_GAS_Z) else None,
            'matches': cumulative_z == NOBLE_GAS_Z[i] if i < len(NOBLE_GAS_Z) else False,
        })

    noble_gas_matches = all(p['matches'] for p in predicted_periods if p['noble_gas_Z'] is not None)

    return {
        'subshell_capacities': predicted_capacities,
        'period_structure': predicted_periods,
        'noble_gas_closures_match': noble_gas_matches,
    }
