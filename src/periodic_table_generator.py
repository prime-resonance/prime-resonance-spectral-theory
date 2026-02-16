"""
Atomic Model Generator — Periodic Table from Primorial Reciprocity Framework.

Generates the full periodic table (Z=1–86) using the 6-step Slater ΔE pipeline,
computes successive ionization energies (IE_1..IE_k), and compares against NIST
reference values.

Output:
    - First ionization energy comparison table (Z=1–86)
    - Successive ionization energy comparison (Z=1–36, up to 10 levels)
    - Error statistics (MAPE, median, max)
    - log(IE_k) vs k data for plotting
"""

import math
from typing import Dict, List, Tuple, Optional

from .atomic_physics.energy import (
    ionization_energy,
    successive_ionization_energies,
)
from .atomic_physics.orbitals import electron_configuration
from .atomic_physics.reporting import config_to_string
from .nist_data import (
    NIST_IONIZATION_ENERGIES,
    get_element_symbol,
)
from .nist_successive_ie import (
    NIST_SUCCESSIVE_IE,
    ELEMENT_SYMBOLS,
    get_nist_successive_ie,
)


def generate_first_ie_table(z_max: int = 86) -> List[Dict]:
    """
    Generate first ionization energy comparison table for Z=1..z_max.
    
    Returns:
        List of dicts with keys: Z, symbol, config, predicted_IE, nist_IE, error_pct
    """
    results = []
    for Z in range(1, z_max + 1):
        symbol = get_element_symbol(Z)
        config = electron_configuration(Z)
        config_str = config_to_string(config)
        predicted = ionization_energy(Z)
        nist_entry = NIST_IONIZATION_ENERGIES.get(Z)
        nist = nist_entry[2] if nist_entry else 0.0
        
        if nist > 0:
            error_pct = 100.0 * (predicted - nist) / nist
        else:
            error_pct = float('nan')
        
        results.append({
            'Z': Z,
            'symbol': symbol,
            'config': config_str,
            'predicted_IE': predicted,
            'nist_IE': nist,
            'error_pct': error_pct,
        })
    
    return results


def compute_ie_statistics(results: List[Dict]) -> Dict:
    """
    Compute error statistics from a first-IE results table.
    """
    errors = [abs(r['error_pct']) for r in results if not math.isnan(r['error_pct'])]
    if not errors:
        return {'mape': 0, 'median': 0, 'max_error': 0, 'count': 0}
    
    errors_sorted = sorted(errors)
    n = len(errors_sorted)
    median = errors_sorted[n // 2] if n % 2 == 1 else (errors_sorted[n // 2 - 1] + errors_sorted[n // 2]) / 2
    
    return {
        'mape': sum(errors) / n,
        'median': median,
        'max_error': max(errors),
        'min_error': min(errors),
        'count': n,
        'within_5pct': sum(1 for e in errors if e <= 5.0),
        'within_10pct': sum(1 for e in errors if e <= 10.0),
        'within_20pct': sum(1 for e in errors if e <= 20.0),
    }


def generate_successive_ie_table(z_max: int = 36, k_max: int = 10) -> List[Dict]:
    """
    Generate successive ionization energy comparison for Z=1..z_max.
    
    Returns:
        List of dicts with keys: Z, symbol, k, predicted_IE, nist_IE, error_pct
    """
    results = []
    for Z in range(1, z_max + 1):
        symbol = ELEMENT_SYMBOLS.get(Z, get_element_symbol(Z))
        nist_ies = get_nist_successive_ie(Z)
        
        k_limit = min(k_max, Z)
        if nist_ies:
            k_limit = min(k_limit, len(nist_ies))
        
        predicted_ies = successive_ionization_energies(Z, max_k=k_limit)
        
        for k in range(len(predicted_ies)):
            nist_val = nist_ies[k] if nist_ies and k < len(nist_ies) else None
            
            if nist_val and nist_val > 0:
                error_pct = 100.0 * (predicted_ies[k] - nist_val) / nist_val
            else:
                error_pct = float('nan')
            
            results.append({
                'Z': Z,
                'symbol': symbol,
                'k': k + 1,
                'predicted_IE': predicted_ies[k],
                'nist_IE': nist_val,
                'error_pct': error_pct,
            })
    
    return results


def compute_successive_ie_statistics(results: List[Dict]) -> Dict:
    """
    Compute error statistics from successive IE results, grouped by ionization level.
    """
    # Overall statistics
    all_errors = [abs(r['error_pct']) for r in results if not math.isnan(r.get('error_pct', float('nan')))]
    
    # Per-level statistics
    by_level = {}
    for r in results:
        k = r['k']
        if k not in by_level:
            by_level[k] = []
        if not math.isnan(r.get('error_pct', float('nan'))):
            by_level[k].append(abs(r['error_pct']))
    
    level_stats = {}
    for k, errors in sorted(by_level.items()):
        if errors:
            n = len(errors)
            level_stats[k] = {
                'mape': sum(errors) / n,
                'median': sorted(errors)[n // 2],
                'max': max(errors),
                'count': n,
            }
    
    overall = {}
    if all_errors:
        n = len(all_errors)
        overall = {
            'mape': sum(all_errors) / n,
            'median': sorted(all_errors)[n // 2],
            'max_error': max(all_errors),
            'count': n,
        }
    
    return {
        'overall': overall,
        'by_level': level_stats,
    }


def generate_log_ie_plot_data(z_list: Optional[List[int]] = None, k_max: int = 10) -> Dict[int, List[Tuple[int, float]]]:
    """
    Generate log(IE_k) vs k data for plotting.
    
    Args:
        z_list: List of Z values to plot (default: selected representative elements)
        k_max: Maximum ionization level
        
    Returns:
        Dict mapping Z -> list of (k, log10(IE_k)) tuples
    """
    if z_list is None:
        z_list = [1, 2, 3, 6, 7, 8, 10, 11, 12, 13, 14, 17, 18, 19, 20, 26, 29, 36]
    
    plot_data = {}
    for Z in z_list:
        ies = successive_ionization_energies(Z, max_k=k_max)
        plot_data[Z] = [(k + 1, math.log10(ie)) for k, ie in enumerate(ies) if ie > 0]
    
    return plot_data


def format_first_ie_markdown(results: List[Dict], stats: Dict) -> str:
    """
    Format first IE results as a Markdown table.
    """
    lines = []
    lines.append("## First Ionization Energy: Framework vs NIST")
    lines.append("")
    lines.append(f"**Statistics (Z=1–{len(results)}):**")
    lines.append(f"- MAPE: {stats['mape']:.2f}%")
    lines.append(f"- Median |error|: {stats['median']:.2f}%")
    lines.append(f"- Max |error|: {stats['max_error']:.2f}%")
    lines.append(f"- Within 5%: {stats['within_5pct']}/{stats['count']}")
    lines.append(f"- Within 10%: {stats['within_10pct']}/{stats['count']}")
    lines.append(f"- Within 20%: {stats['within_20pct']}/{stats['count']}")
    lines.append("")
    lines.append("| Z | Symbol | Config | Predicted (eV) | NIST (eV) | Error (%) |")
    lines.append("|---|--------|--------|---------------:|----------:|----------:|")
    
    for r in results:
        err_str = f"{r['error_pct']:+.2f}" if not math.isnan(r['error_pct']) else "N/A"
        lines.append(
            f"| {r['Z']:>2} | {r['symbol']:>4} | {r['config']:<20} | "
            f"{r['predicted_IE']:>14.3f} | {r['nist_IE']:>9.3f} | {err_str:>9} |"
        )
    
    return "\n".join(lines)


def format_successive_ie_markdown(results: List[Dict], stats: Dict) -> str:
    """
    Format successive IE results as a Markdown table.
    """
    lines = []
    lines.append("## Successive Ionization Energies: Framework vs NIST")
    lines.append("")
    
    overall = stats.get('overall', {})
    if overall:
        lines.append(f"**Overall Statistics:**")
        lines.append(f"- MAPE: {overall.get('mape', 0):.2f}%")
        lines.append(f"- Median |error|: {overall.get('median', 0):.2f}%")
        lines.append(f"- Total comparisons: {overall.get('count', 0)}")
        lines.append("")
    
    lines.append("**Per-Level MAPE:**")
    lines.append("")
    lines.append("| Level k | MAPE (%) | Median (%) | Max (%) | Count |")
    lines.append("|--------:|---------:|-----------:|--------:|------:|")
    for k, ls in sorted(stats.get('by_level', {}).items()):
        lines.append(f"| {k:>7} | {ls['mape']:>8.2f} | {ls['median']:>10.2f} | {ls['max']:>7.1f} | {ls['count']:>5} |")
    lines.append("")
    
    lines.append("| Z | Symbol | k | Predicted (eV) | NIST (eV) | Error (%) |")
    lines.append("|---|--------|---|---------------:|----------:|----------:|")
    
    for r in results:
        nist_str = f"{r['nist_IE']:.3f}" if r['nist_IE'] is not None else "—"
        err_str = f"{r['error_pct']:+.2f}" if not math.isnan(r.get('error_pct', float('nan'))) else "—"
        lines.append(
            f"| {r['Z']:>2} | {r['symbol']:>4} | {r['k']} | "
            f"{r['predicted_IE']:>14.3f} | {nist_str:>9} | {err_str:>9} |"
        )
    
    return "\n".join(lines)


def run_full_comparison(z_max_first: int = 86, z_max_successive: int = 36,
                        k_max: int = 10) -> Dict:
    """
    Run the full periodic table generation and NIST comparison.
    
    Returns dict with all results and formatted output.
    """
    # First IE
    first_ie_results = generate_first_ie_table(z_max_first)
    first_ie_stats = compute_ie_statistics(first_ie_results)
    
    # Successive IE
    successive_ie_results = generate_successive_ie_table(z_max_successive, k_max)
    successive_ie_stats = compute_successive_ie_statistics(successive_ie_results)
    
    # log(IE) plot data
    plot_data = generate_log_ie_plot_data(k_max=k_max)
    
    # Format markdown
    first_md = format_first_ie_markdown(first_ie_results, first_ie_stats)
    successive_md = format_successive_ie_markdown(successive_ie_results, successive_ie_stats)
    
    return {
        'first_ie': {
            'results': first_ie_results,
            'stats': first_ie_stats,
            'markdown': first_md,
        },
        'successive_ie': {
            'results': successive_ie_results,
            'stats': successive_ie_stats,
            'markdown': successive_md,
        },
        'plot_data': plot_data,
    }


if __name__ == '__main__':
    print("=" * 80)
    print("PRIMORIAL RECIPROCITY FRAMEWORK — PERIODIC TABLE GENERATOR")
    print("Generating periodic table and comparing against NIST values...")
    print("=" * 80)
    print()
    
    output = run_full_comparison()
    
    print(output['first_ie']['markdown'])
    print()
    print(output['successive_ie']['markdown'])
