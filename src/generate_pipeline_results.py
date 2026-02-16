
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Add src to python path so imports work
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src.atomic_physics.energy import successive_ionization_energies
from src.nist_successive_ie import NIST_SUCCESSIVE_IE, ELEMENT_SYMBOLS
from src.nist_data import NIST_IONIZATION_ENERGIES

def get_experimental_ie(Z, k):
    """
    Retrieve experimental ionization energy for element Z at level k.
    k is 1-indexed (1st IE, 2nd IE, etc.)
    """
    # Check NIST_SUCCESSIVE_IE first
    if Z in NIST_SUCCESSIVE_IE:
        ies = NIST_SUCCESSIVE_IE[Z]
        if k <= len(ies):
            return ies[k-1]
            
    # Fallback to NIST_IONIZATION_ENERGIES for k=1
    if k == 1 and Z in NIST_IONIZATION_ENERGIES:
        return NIST_IONIZATION_ENERGIES[Z][2]
        
    return None

def main():
    results = []
    
    # Range of Z to process
    Z_min = 1
    Z_max = 86
    
    output_dir = "output"
    plots_dir = os.path.join(output_dir, "plots")
    
    print(f"Running pipeline for Z={Z_min} to {Z_max}...")
    
    for Z in range(Z_min, Z_max + 1):
        max_k = min(10, Z)
        
        # Calculate theoretical values
        try:
            theoretical_ies = successive_ionization_energies(Z, max_k=max_k)
        except Exception as e:
            print(f"Error calculating for Z={Z}: {e}")
            theoretical_ies = []
            
        # Collect data for this Z
        ks = []
        theo_vals = []
        exp_vals = []
        
        for k in range(1, max_k + 1):
            if k - 1 < len(theoretical_ies):
                theo = theoretical_ies[k-1]
            else:
                theo = None # Should not happen if calculation successful
            
            exp = get_experimental_ie(Z, k)
            
            # Record result
            record = {
                'Z': Z,
                'Symbol': ELEMENT_SYMBOLS.get(Z, NIST_IONIZATION_ENERGIES.get(Z, ('?',))[0]),
                'k': k,
                'Experimental_IE': exp,
                'Theoretical_IE': theo
            }
            
            if exp is not None and theo is not None:
                error = theo - exp
                pct_error = (error / exp) * 100
                record['Error_eV'] = error
                record['Error_Pct'] = pct_error
                exp_vals.append(exp)
            else:
                record['Error_eV'] = None
                record['Error_Pct'] = None
                exp_vals.append(np.nan) # Keep placeholder for plotting
                
            results.append(record)
            
            if theo is not None:
                theo_vals.append(theo)
                ks.append(k)
        
        # Generate Plot for this Z
        if theo_vals:
            plt.figure(figsize=(8, 6))
            plt.plot(ks, np.log10(theo_vals), 'bo-', label='Theoretical')
            
            # Plot experimental if we have valid points
            valid_exp_indices = [i for i, v in enumerate(exp_vals) if not np.isnan(v)]
            if valid_exp_indices:
                valid_ks = [ks[i] for i in valid_exp_indices]
                valid_exps = [exp_vals[i] for i in valid_exp_indices]
                plt.plot(valid_ks, np.log10(valid_exps), 'rx--', label='Experimental (NIST)')
                
            plt.xlabel('Ionization Level k')
            plt.ylabel('log10(IE) [eV]')
            plt.title(f'Successive Ionization Energies for Z={Z} ({record["Symbol"]})')
            plt.legend()
            plt.grid(True)
            
            plot_path = os.path.join(plots_dir, f"Z{Z:03d}_IE_plot.png")
            plt.savefig(plot_path)
            plt.close()
            
    # Save Master Error Table
    df = pd.DataFrame(results)
    csv_path = os.path.join(output_dir, "master_error_table.csv")
    df.to_csv(csv_path, index=False)
    print(f"Master error table saved to {csv_path}")
    print(f"Plots saved to {plots_dir}")

if __name__ == "__main__":
    main()
