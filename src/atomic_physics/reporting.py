"""
Reporting Functions for Atomic Model.

Generates periodic table comparisons, accuracy statistics,
alkali metal trends, and orbital structure validation.
"""

from typing import List, Dict, Any, Optional
from ..nist_data import (
    get_ionization_energy, get_element_symbol, get_period,
    NOBLE_GAS_Z, ALKALI_Z, PERIOD_LENGTHS, NIST_IONIZATION_ENERGIES,
    SHELL_FILLING_ORDER
)
from .orbitals import electron_configuration, validate_orbital_structure
from .energy import ionization_energy


def config_to_string(config: List) -> str:
    """Convert electron configuration to string notation."""
    l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    parts = []
    for n, l, count in config:
        parts.append(f"{n}{l_names.get(l, '?')}{count}")
    return ' '.join(parts)


def generate_periodic_table(params: Dict[str, Any] = None) -> List[Dict]:
    """
    Generate periodic table data with predicted vs NIST ionization energies.
    
    Returns list of dicts: Z, symbol, config, predicted_IE, nist_IE, error_pct
    """
    results = []
    
    for Z in range(1, 119):
        nist_ie = get_ionization_energy(Z)
        if nist_ie is None:
            continue
            
        symbol = get_element_symbol(Z)
        config = electron_configuration(Z)
        config_str = config_to_string(config)
        pred_ie = ionization_energy(Z, params)
        
        error_pct = abs(pred_ie - nist_ie) / nist_ie * 100 if nist_ie > 0 else 0.0
        
        results.append({
            'Z': Z,
            'symbol': symbol,
            'config': config_str,
            'predicted_IE': pred_ie,
            'nist_IE': nist_ie,
            'error_pct': error_pct
        })
        
    return results


def generate_periodic_table_summary(params: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Generate summary statistics for the periodic table predictions.
    
    Returns overall and per-period accuracy statistics.
    """
    table = generate_periodic_table(params)
    
    # Overall stats
    errors = [row['error_pct'] for row in table]
    
    summary = {
        'overall': {
            'count': len(errors),
            'mean_error': sum(errors) / len(errors) if errors else 0,
            'max_error': max(errors) if errors else 0,
            'min_error': min(errors) if errors else 0,
            'median_error': sorted(errors)[len(errors) // 2] if errors else 0,
        },
        'periods': {},
        'ranges': {}
    }
    
    # Per-period stats
    for period_num in range(1, 8):
        period_errors = []
        for row in table:
            if get_period(row['Z']) == period_num:
                period_errors.append(row['error_pct'])
        if period_errors:
            summary['periods'][period_num] = {
                'count': len(period_errors),
                'mean_error': sum(period_errors) / len(period_errors),
                'max_error': max(period_errors),
            }
    
    # Range stats (Z=1..36, Z=1..86)
    for label, z_max in [('Z1_36', 36), ('Z1_86', 86)]:
        range_errors = [r['error_pct'] for r in table if r['Z'] <= z_max]
        if range_errors:
            summary['ranges'][label] = {
                'count': len(range_errors),
                'mean_error': sum(range_errors) / len(range_errors),
                'max_error': max(range_errors),
            }
    
    return summary


def alkali_ionization_trend(params: Dict[str, Any] = None) -> List[Dict]:
    """
    Compare predicted vs NIST ionization energies for alkali metals.
    
    IE should decrease monotonically: H > Li > Na > K > Rb > Cs
    """
    results = []
    for Z in ALKALI_Z:
        nist_ie = get_ionization_energy(Z)
        if nist_ie is None:
            continue
        symbol = get_element_symbol(Z)
        pred_ie = ionization_energy(Z, params)
        error_pct = abs(pred_ie - nist_ie) / nist_ie * 100 if nist_ie > 0 else 0.0
        
        results.append({
            'Z': Z,
            'symbol': symbol,
            'predicted_IE': pred_ie,
            'nist_IE': nist_ie,
            'error_pct': error_pct
        })
        
    return results
