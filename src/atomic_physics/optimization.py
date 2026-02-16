"""
Grid Search Optimization for Atomic Model Parameters.

Uses scipy.optimize.differential_evolution for global optimization
of the Slater ΔE pipeline parameters, minimizing Mean Absolute
Percentage Error (MAPE) against NIST data for Z=1..86.

The objective penalizes both average error AND large outliers to
prevent the optimizer from sacrificing a few elements for overall gain.
"""

import numpy as np
from scipy.optimize import differential_evolution
from typing import Dict, Any, List, Tuple
from ..nist_data import get_ionization_energy
from .energy import ionization_energy
from .constants import (
    N_STAR_VALUES, SHIELDING_COEFFS,
    C_PAIR, C_EXCH, C_REL_SP, C_REL_DF
)


def unpack_parameters(x: np.ndarray) -> Dict[str, Any]:
    """
    Convert flat parameter vector to structured params dict.

    Layout (31 params):
      [0..6]   n* for s-orbitals: n=1..7
      [7..11]  n* for p-orbitals: n=2..6
      [12..15] n* for d-orbitals: n=3..6
      [16..18] n* for f-orbitals: n=4..6
      [19..21] s shielding: same_group, one_below, deep_core
      [22..23] p shielding: same_group, one_below
      [24]     df same_group shielding
      [25..26] df_source_reduction: for s, for p
      [27..30] C_PAIR, C_EXCH, C_REL_SP, C_REL_DF
    """
    n_star = {
        1: {0: x[0]},
        2: {0: x[1], 1: x[7]},
        3: {0: x[2], 1: x[8], 2: x[12]},
        4: {0: x[3], 1: x[9], 2: x[13], 3: x[16]},
        5: {0: x[4], 1: x[10], 2: x[14], 3: x[17]},
        6: {0: x[5], 1: x[11], 2: x[15], 3: x[18]},
        7: {0: x[6], 1: x[6], 2: x[6], 3: x[6]},
    }

    shielding = {
        's': {
            'same_group': x[19],
            'one_below': x[20],
            'deep_core': x[21],
        },
        'p': {
            'same_group': x[22],
            'one_below': x[23],
            'deep_core': 1.00,
        },
        'df': {
            'same_group': x[24],
            'one_below': 1.00,
            'deep_core': 1.00,
        },
        'df_source_reduction': {
            's': x[25],
            'p': x[26],
        },
    }

    return {
        'n_star': n_star,
        'shielding_coeffs': shielding,
        'c_pair': x[27],
        'c_exch': x[28],
        'c_rel_sp': x[29],
        'c_rel_df': x[30],
    }


def pack_parameters() -> Tuple[np.ndarray, List[Tuple[float, float]]]:
    """
    Pack current default parameters into flat vector and bounds.
    """
    x0 = np.array([
        # n* s [0..6]
        1.0, 2.0, 3.0, 3.50, 4.00, 4.15, 4.30,
        # n* p [7..11]
        2.0, 3.0, 3.55, 4.00, 4.20,
        # n* d [12..15]
        2.80, 3.45, 4.00, 4.20,
        # n* f [16..18]
        3.40, 3.80, 4.00,
        # shielding [19..26]
        0.30, 0.85, 0.95, 0.35, 0.93, 0.36, 0.60, 0.75,
        # corrections [27..30]
        0.45, 0.35, 0.08, 0.02,
    ])

    bounds = [
        # n* s
        (0.90, 1.10), (1.7, 2.3), (2.4, 3.3), (3.0, 4.0),
        (3.3, 4.5), (3.5, 5.0), (3.8, 5.2),
        # n* p
        (1.7, 2.4), (2.4, 3.3), (3.0, 4.0), (3.3, 4.5), (3.5, 5.0),
        # n* d
        (2.0, 3.5), (2.8, 4.0), (3.0, 4.5), (3.3, 5.0),
        # n* f
        (2.8, 4.2), (3.2, 4.5), (3.5, 5.0),
        # shielding
        (0.20, 0.40), (0.60, 1.00), (0.80, 1.05),
        (0.25, 0.45), (0.70, 1.00),
        (0.25, 0.50),
        (0.30, 1.00), (0.40, 1.00),
        # corrections
        (0.0, 2.0), (0.0, 2.0), (0.0, 0.3), (0.0, 0.15),
    ]

    return x0, bounds


def objective_function(x: np.ndarray) -> float:
    """
    Composite objective: MAPE + penalty for large outliers.

    obj = mean(errors) + 0.3 * max(errors)

    This prevents the optimizer from sacrificing individual elements
    (like Sc, Ti) to get slightly better average performance elsewhere.
    """
    params = unpack_parameters(x)

    errors = []
    for Z in range(1, 87):
        nist_val = get_ionization_energy(Z)
        if nist_val is None or nist_val <= 0:
            continue

        try:
            pred_val = ionization_energy(Z, params)
        except Exception:
            errors.append(200.0)
            continue

        if pred_val <= 0.01:
            errors.append(200.0)
            continue

        error_pct = abs(pred_val - nist_val) / nist_val * 100
        errors.append(error_pct)

    if not errors:
        return 1000.0

    mean_err = np.mean(errors)
    max_err = np.max(errors)

    # Composite: penalize both mean and worst-case
    return mean_err + 0.3 * max_err


def optimize_parameters(maxiter: int = 100, seed: int = 42,
                        popsize: int = 30, tol: float = 0.01) -> Tuple[Dict[str, Any], float]:
    """
    Run differential evolution to find optimal parameters.

    Returns:
        (best_params_dict, best_mape)
    """
    x0, bounds = pack_parameters()

    print(f"Starting differential evolution optimization...")
    print(f"  Parameters: {len(x0)}, Population: {popsize}, MaxIter: {maxiter}")

    result = differential_evolution(
        objective_function,
        bounds,
        x0=x0,
        maxiter=maxiter,
        popsize=popsize,
        seed=seed,
        tol=tol,
        mutation=(0.5, 1.0),
        recombination=0.7,
        disp=True,
        polish=True,
    )

    print(f"\nOptimization finished. Success: {result.success}")
    print(f"Composite objective: {result.fun:.4f}")

    best_params = unpack_parameters(result.x)

    # Calculate actual MAPE
    errors = []
    for Z in range(1, 87):
        nist_val = get_ionization_energy(Z)
        if nist_val is None or nist_val <= 0:
            continue
        pred_val = ionization_energy(Z, best_params)
        error_pct = abs(pred_val - nist_val) / nist_val * 100
        errors.append(error_pct)

    mape = np.mean(errors) if errors else 0
    print(f"Actual MAPE: {mape:.4f}%")
    print(f"Max error: {np.max(errors):.4f}%")

    # Print optimized values
    print("\n--- Optimized n* Values ---")
    ns = best_params['n_star']
    l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    for n in sorted(ns.keys()):
        vals = ', '.join(f"{l_names[l]}={v:.4f}" for l, v in sorted(ns[n].items()))
        print(f"  n={n}: {vals}")

    sc = best_params['shielding_coeffs']
    print("\n--- Optimized Shielding ---")
    for key in ['s', 'p', 'df']:
        print(f"  {key}: {sc[key]}")
    if 'df_source_reduction' in sc:
        print(f"  df_source_reduction: {sc['df_source_reduction']}")

    print(f"\n--- Correction Coefficients ---")
    print(f"  C_PAIR:   {best_params['c_pair']:.4f}")
    print(f"  C_EXCH:   {best_params['c_exch']:.4f}")
    print(f"  C_REL_SP: {best_params['c_rel_sp']:.4f}")
    print(f"  C_REL_DF: {best_params['c_rel_df']:.4f}")

    # Print raw vector for embedding
    print(f"\n--- Raw x vector ---")
    print(f"x = {result.x.tolist()}")

    return best_params, mape


if __name__ == '__main__':
    optimize_parameters()
