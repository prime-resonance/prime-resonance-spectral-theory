# Primorial Reciprocity Framework

**Deriving Standard Model Constants and Atomic Scales from the Arithmetic of 2310 = 2 × 3 × 5 × 7 × 11**

[![Tests](https://img.shields.io/badge/tests-383%20passing-brightgreen)]()
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

---

## Overview

This repository contains a computationally verified framework demonstrating that all dimensionless mass ratios, coupling constants, mixing angles of the Standard Model, and atomic ionization energies can be expressed through one structural principle: the decomposition of the primorial **2310 = 2 × 3 × 5 × 7 × 11** into reciprocity channels (quadratic, cubic, quintic, plus higher extensions).

The generalized twist unit **T(P₅) = 6480 = 2⁴ × 3⁴ × 5** serves as the core scale, with the sub-primorial identity T(P₃) = 108 = 2² × 3³ as the mass quantum.

### Key Results: Particle Physics

| Constant | Formula | Predicted | Measured | Error |
|----------|---------|-----------|----------|-------|
| m_p/m_e | 17 × 108 | 1836 | 1836.15 | 0.008% |
| α⁻¹ | 108 + 29 + 1/27 | 137.037 | 137.036 | 0.0007% |
| m_H/m_e | 21×108² − 27 | 244917 | 245107 | 0.08% |
| m_τ/m_e | 31×108 + 125 + 4 | 3477 | 3477.2 | 0.007% |
| **Δm²₃₂** | from 108³×8×totient | **2.35e-3** | **2.453e-3** | **4.2%** |
| sin²θ_W | 3/13 | 0.2308 | 0.2312 | 0.2% |
| sin θ_C | 29/128 | 0.2266 | 0.2257 | 0.3% |
| sin²θ_12 | 4/13 | 0.3077 | 0.304 | 1.2% |
| sin²θ_23 | 4/7 | 0.5714 | 0.570 | 0.2% |
| sin²θ_13 | 1/45 | 0.0222 | 0.0222 | 0.1% |

### Key Results: Atomic Physics

The same primorial channels govern atomic ionization energies via a 6-step Slater ΔE pipeline:

| Metric | First IE (Z=1–86) | Successive IE (Z=1–36) |
|--------|-------------------|------------------------|
| MAPE | **5.74%** | **20.59%** |
| Median | 4.26% | 10.78% |
| Within 5% | 52/86 (60%) | — |
| Within 10% | 68/86 (79%) | — |
| Comparisons | 86 | 315 |

Selected predictions (first IE):

| Element | Z | Predicted (eV) | NIST (eV) | Error |
|---------|---|---------------|-----------|-------|
| F | 9 | 17.42 | 17.42 | 0.04% |
| Y | 39 | 6.21 | 6.22 | 0.07% |
| I | 53 | 10.44 | 10.45 | 0.14% |
| Rn | 86 | 10.69 | 10.75 | 0.54% |
| N | 7 | 14.43 | 14.53 | 0.71% |
| Ni | 28 | 7.58 | 7.64 | 0.83% |
| Lanthanides (Ce–Yb) | 58–70 | — | — | **2.4% mean** |

See the [full particle physics paper](papers/primorial-reciprocity-and-the-mass-spectrum.md) and [atomic physics paper](papers/atomic-ionization-from-primorial-reciprocity.md) for all derivations.

## Note on Units and the Higgs Mass

While the framework primarily predicts **dimensionless ratios** (masses relative to the electron), the Higgs mass value of ~125 GeV numerically resonates with the quintic channel power ($5^3 = 125$). This leads to a dual consistency:
1.  **Dimensionless**: $m_H/m_e \approx 21 \times 108^2 - 27 \approx 244917$ (matches experiment within 0.08%).
2.  **Scale**: The electroweak scale appears to align such that the quintic channel's "natural unit" coincides with the GeV.

## Quick Start

```bash
# Clone
git clone https://github.com/sschepis/prime-resonance-spectral-theory.git
cd prime-resonance-spectral-theory

# Install dependencies
pip install numpy scipy pytest

# Run all 383 tests
pytest tests/ -v

# Generate the periodic table with NIST comparison
python3.10 -c "
from src.periodic_table_generator import run_full_comparison
output = run_full_comparison()
print(output['first_ie']['markdown'])
print(output['successive_ie']['markdown'])
"
```

## The Framework

### Central Objects

| Object | Definition | Role |
|--------|-----------|------|
| 2310 = 2×3×5×7×11 | Full primorial | Sieve modulus |
| 6480 = 2⁴×3⁴×5 | Twist unit T(P₅) | Mass quantum (generalized) |
| 108 = 2²×3³ | Sub-twist T(P₃) | Base mass quantum |
| 125 = 5³ | Quintic power | Higgs mass (GeV) |
| 8 = φ(30) | Euler totient (base) | Coprime multiplicity |
| 480 = φ(2310) | Euler totient (full) | Full coprime multiplicity |

### The Three-Channel Framework

The primorial 30 = 2 × 3 × 5 produces three independent reciprocity channels:

| Channel | Prime | Reciprocity | Ring | Contribution | Physical Role |
|---------|-------|-------------|------|-------------|---------------|
| Quadratic | 2 | Legendre | Z | 2² = 4 | Sub-shell structure, non-Hermiticity |
| Cubic | 3 | Eisenstein | Z[ω] | 3³ = 27 | Mass corrections {±3^k}, 3-adic tower |
| Quintic | 5 | Cyclotomic | Z[ζ₅] | 5³ = 125 | Higgs boson mass (GeV) |

**Combined:** 108 = 2² × 3³ is the mass quantum. Particle masses follow:
- **Positive branch:** m/m_e = n × 108 ± 3^k (proton, muon, tau, W, Z, t, b, c, Higgs)
- **Inverse branch:** m_ν = m_e / (108³ × 8 × 3^k) (neutrinos)

Adding primes 7 and 11 multiplies the twist unit to **6480** and refines every channel.

### Channel Mapping to Atomic Orbitals

| *l* | Label | Prime | Reciprocity Channel | Capacity | Atomic Role |
|-----|-------|-------|---------------------|----------|-------------|
| 0 | s | 2 | Quadratic | 2 | Penetrating, tightest binding |
| 1 | p | 3 | Cubic | 6 | Intermediate, pairing effects |
| 2 | d | 5 | Quintic | 10 | Compact, exchange stabilization |
| 3 | f | 7 | Septic | 14 | Most compact, lanthanide accuracy |

### Mixing Angles

The framework derives mixing angles from ratios of small primes and twist-related numbers:

*   **Cabibbo Angle**: $\sin \theta_C \approx 29/128$ (29 is the boundary prime, $128=2^7$)
*   **Weinberg Angle**: $\sin^2 \theta_W \approx 3/13$
*   **Neutrino Mixing (PMNS)**:
    *   $\sin^2 \theta_{12} \approx 4/13$
    *   $\sin^2 \theta_{23} \approx 4/7$
    *   $\sin^2 \theta_{13} \approx 1/45$

## The Atomic Model Pipeline

A 6-step Slater ΔE pipeline predicts ionization energies with 31 globally-optimized parameters and a parameter-free 3-adic tower correction:

1. **Electron Configuration** — Aufbau filling order
2. **Effective Quantum Numbers** n*(n,l) — 19 optimized values encoding orbital penetration
3. **l-Dependent Shielding** — 9 coefficients + 2 source-reduction factors + 3-adic tower
4. **Pairing Correction** — Beyond-half-filled subshell repulsion (C_pair = 0.139)
5. **Exchange Stabilization** — Half-filled subshell stability (C_exch = 1.457)
6. **Relativistic Corrections** (Z > 36) — s/p contraction, d/f expansion

The **3-adic tower correction** δ = −1/3^(n−2) for valence s-orbitals (n > 2) introduces zero new free parameters and reduces Na error from 39.3% → 4.8%.

## Repository Structure

```
prime-resonance-spectral-theory/
├── src/                              # Source modules
│   ├── prime_utils.py                # Sieve, Legendre/Kronecker symbols, 2310 constants
│   ├── mod30_structure.py            # Period-24 cycles, residue classification
│   ├── legendre_network.py           # Interaction matrices, qubit scoring
│   ├── spectral_analysis.py          # Eigenspectra, modularity, random comparison
│   ├── resonance_annealer.py         # Topological qubits, CSP solver
│   ├── periodicity.py                # FFT periodicity, 108 analysis
│   ├── collapse_hamiltonian.py       # Non-Hermitian Hamiltonian, collapse dynamics
│   ├── feigenbaum_analysis.py        # Spectral determinant, bifurcation
│   ├── spectral_constants.py         # Scaling laws, Wigner surmise, traces
│   ├── analytic_108.py               # CRT proof, General Twist Formula
│   ├── primorial_spectrum.py         # Inter-primorial analysis, mass search
│   ├── cubic_residue.py              # Cubic reciprocity, Eisenstein integers
│   ├── quintic_residue.py            # Quintic channel, Higgs derivation
│   ├── neutrino_prediction.py        # Inverse twist tower, Δm²₃₂ prediction
│   ├── mixing_angles.py              # P5 scale physics, fine structure
│   ├── nist_data.py                  # NIST first IE data for Z=1–118
│   ├── nist_successive_ie.py         # NIST successive IE data for Z=1–36
│   ├── atomic_model.py               # Facade for atomic physics pipeline
│   ├── periodic_table_generator.py   # Full periodic table generator + NIST comparison
│   └── atomic_physics/               # Modular atomic physics subpackage
│       ├── constants.py              # 31 optimized parameters
│       ├── orbitals.py               # Electron configuration, quantum numbers
│       ├── shielding.py              # l-dependent shielding, 3-adic tower
│       ├── energy.py                 # Total energy, first + successive IE
│       ├── optimization.py           # Differential evolution optimizer
│       └── reporting.py              # Periodic table generation, statistics
│
├── tests/                            # Test suite (383 tests)
│   ├── test_mod30_structure.py
│   ├── test_legendre_asymmetry.py
│   ├── test_spectral_properties.py
│   ├── test_108_periodicity.py
│   ├── test_resonance_annealing.py
│   ├── test_digital_root_cycles.py
│   ├── test_collapse_dynamics.py
│   ├── test_spectral_constants.py
│   ├── test_analytic_and_predictions.py
│   ├── test_primorial_spectrum.py
│   ├── test_cubic_residue.py
│   ├── test_quintic_residue.py
│   ├── test_neutrino_prediction.py
│   ├── test_mixing_angles.py
│   ├── test_2310_framework.py        # 58 tests for 2310 primorial
│   ├── test_atomic_model.py          # 32 tests for atomic model
│   └── test_periodic_table_generator.py  # 47 tests for periodic table + successive IE
│
├── papers/
│   ├── primorial-reciprocity-and-the-mass-spectrum.md
│   ├── legendre-weighted-prime-hamiltonian.md
│   ├── general-twist-formula-and-mass-quantization.md
│   └── atomic-ionization-from-primorial-reciprocity.md
├── pyproject.toml
└── README.md
```

## Key Theorems (Proven)

### Theorem 1: The 108 Identity
The sum of one period-24 digital root cycle of the coprime-to-30 sequence is exactly 108 = 2² × 3³.

### Theorem 2: The 6480 Identity
The sum of one full cycle of digital roots of numbers coprime to 2310 is exactly **6480 = 2⁴ × 3⁴ × 5**. This is T(P₅) from the General Twist Formula.

### Theorem 3: The General Twist Formula
T(P_k) = 3³ × Π_{p|P_k, p≥5} (p−1). Consecutive ratios are Euler totients: φ(5)=4, φ(7)=6, φ(11)=10.

### Theorem 4: Non-Hermiticity from Quadratic Reciprocity
The Legendre-weighted Hamiltonian is intrinsically non-Hermitian because (p/q) ≠ (q/p) for primes p,q ≡ 3 (mod 4).

### Theorem 5: The 288 Connection
288 = T(P₃) × φ(P₃) / 3, connecting the figure-eight knot's PSL(2,Z[ω]) symmetry to the twist unit.

## Falsifiable Predictions

### Particle Physics
1. **Neutrino mass ratio m₃/m₂ = 3 exactly** — testable at JUNO, DUNE, Hyper-Kamiokande
2. **Σm_ν ≈ 72.8 meV** — testable via CMB-S4 and galaxy surveys
3. **Δm²₃₂ ≈ 2.35 × 10⁻³ eV²** — currently measured at 2.453 × 10⁻³ (4.2% difference)
4. **P₆ = 30030 scale** predicts next B-meson family or new resonances

### Atomic Physics
5. **No intermediate states violating the n × 108 ± 3^k × totient-correction rule**
6. **3-adic tower correction** −1/3^(n−2) improves all alkali s-orbital predictions
7. **Lanthanide IE accuracy < 3% MAPE** with no per-element fitting

## Test Suite Summary

| Test Suite | Tests | Coverage |
|-----------|-------|----------|
| Original PRF tests | 246 | Mass spectrum, twist units, spectral analysis |
| 2310 framework tests | 58 | 6480 identity, P₅ physics, refined predictions |
| Atomic model tests | 32 | Orbital structure, shielding, IE accuracy, trends |
| Periodic table generator | 47 | Successive IE, NIST comparison, statistics |
| **Total** | **383** | **All passing** |

## Citation

```bibtex
@article{schepis2026primorial,
  title={Primorial Reciprocity Framework: Deriving Standard Model
         Constants and Atomic Scales from the Arithmetic of 2310},
  author={Schepis, Sebastian},
  year={2026},
  note={383 computational tests. Code: github.com/sschepis/prime-resonance-spectral-theory}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
