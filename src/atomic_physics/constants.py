"""
Atomic Physics Constants and Parameters.

All tunable parameters for the Justified Slater ΔE Pipeline.
Optimized by differential evolution (scipy) to minimize MAPE
against NIST ionization energies for Z=1..86.

After optimization: MAPE ≈ 8.8% for Z=1..86, ~14% for Z=1..36.
"""

# ─── Fundamental Constants ────────────────────────────────────────────
RYDBERG_EV = 13.605693  # eV (NIST CODATA 2018)

# ─── Step 2: Effective Quantum Numbers n*(n, l) ───────────────────────
# Optimized via differential evolution. Physical justification:
# - s-orbitals penetrate core → lower n* → tighter binding
# - d,f orbitals are compact → lower n* than standard Slater
# - 6s has relativistic contraction → reduced n*
# - p orbitals have less penetration than s → slightly different n*
N_STAR_VALUES = {
    1: {0: 0.9640},           # 1s: near 1.0, slight optimization
    2: {0: 1.9019, 1: 1.9134},  # 2s/2p: near 2.0
    3: {0: 2.7608, 1: 2.5668, 2: 2.7057},  # 3s < 3.0 (penetration), 3p < 3s
    4: {0: 3.0646, 1: 3.4237, 2: 3.1125, 3: 3.8234},  # 4s most penetrating
    5: {0: 3.5961, 1: 3.7130, 2: 3.6245, 3: 4.0717},
    6: {0: 4.6594, 1: 4.1501, 2: 4.2449, 3: 4.5365},  # 6s: relativistic
    7: {0: 4.9058, 1: 4.9058, 2: 4.9058, 3: 4.9058},
}

# ─── Step 3: l-Dependent Shielding Constants ─────────────────────────
# Optimized coefficients for shielding contributions.
#
# Physical justification:
# - s penetrates core → inner electrons shield s LESS (one_below ~0.97)
# - p doesn't penetrate as deep → shielded more (one_below ~0.98)
# - d/f are compact → ALL lower groups shield them fully (~1.00)
# - same_group coefficients reflect mutual shielding within a subshell
# - df_source_reduction: d/f source orbitals are compact, so they
#   shield outer s/p less effectively (factor ~0.94-0.99)
SHIELDING_COEFFS = {
    's': {
        'same_group': 0.2852,   # Same (ns) group: weak mutual shielding
        'one_below': 0.9660,    # (n-1) shell: optimized
        'deep_core': 0.9748,    # (n-2) and below: near-full shielding
    },
    'p': {
        'same_group': 0.3393,   # Same (np) group
        'one_below': 0.9757,    # (n-1) shell: p shielded more than s
        'deep_core': 1.00,      # Deep core: full shielding for p
    },
    'df': {
        'same_group': 0.4708,   # Same (nd/nf) group: higher mutual shielding
        'one_below': 1.00,      # Inner shells fully shield d/f
        'deep_core': 1.00,      # Deep core: full shielding
    },
    # d/f source orbitals are compact, so they shield s/p targets less.
    'df_source_reduction': {
        's': 0.9875,   # d shields s slightly less
        'p': 0.9430,   # d shields p moderately less
    },
}

# ─── Step 4: Pairing Correction ──────────────────────────────────────
# For l >= 1 subshells beyond half-filling, paired electrons experience
# Coulomb repulsion → easier ionization (lower IE).
# Δ_pair = -C_PAIR × (n_e - half) / half  [eV]
# Optimized value is smaller than initial (0.45) because the
# improved shielding already captures some pairing effects.
C_PAIR = 0.1390

# ─── Step 5: Exchange Stabilization ──────────────────────────────────
# Exactly half-filled subshells (p³, d⁵, f⁷) have maximum exchange energy.
# Removing one electron costs extra energy.
# Δ_exchange = +C_EXCH × (l / 2)  [eV]
# Optimized to be larger than initial (0.35) to better capture
# the N > O and Mn stability effects.
C_EXCH = 1.4567

# ─── Step 6: Relativistic Corrections (Z > 36) ──────────────────────
# Using α = 1/137.036:
# s/p contraction: IE *= (1 + C_REL_SP × Z²α²)
# d/f expansion:   IE *= (1 - C_REL_DF × Z²α²)
C_REL_SP = 0.1038
C_REL_DF = 0.0513
