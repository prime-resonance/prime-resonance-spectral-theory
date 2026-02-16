"""
Cubic Residue Symbol and the 3-Adic Correction Tower.

The Legendre (quadratic) symbol governs the 2-adic (mod-4) structure:
    (p/q)(q/p) = (-1)^{(p-1)/2·(q-1)/2}  → non-Hermiticity

The CUBIC residue symbol should govern the 3-adic (mod-9) structure:
    χ₃(a, p) = a^{(p-1)/3} mod p  for p ≡ 1 (mod 3)

This module implements the cubic residue symbol over the Eisenstein integers
Z[ω] where ω = e^{2πi/3} = (-1 + i√3)/2, and tests whether it generates
the {0, ±3, ±9, ±27} correction tower observed in particle masses.

KEY HYPOTHESIS: Just as Legendre asymmetry (2-adic) → non-Hermiticity,
cubic residue asymmetry (3-adic) → the correction hierarchy.
The COMBINED (2×3)-adic structure gives 108 = 2²×3³ as the fundamental unit.
"""

from typing import Dict, List, Tuple, Optional
import math
import numpy as np


# ============================================================================
# Eisenstein Integers Z[ω]
# ============================================================================

# ω = e^{2πi/3} = (-1 + i√3)/2
OMEGA_REAL = -0.5
OMEGA_IMAG = math.sqrt(3) / 2


class EisensteinInt:
    """
    An Eisenstein integer a + bω where ω = e^{2πi/3}.

    Norm: N(a + bω) = a² - ab + b²
    """

    def __init__(self, a: int, b: int):
        self.a = a
        self.b = b

    @property
    def norm(self) -> int:
        """Norm: a² - ab + b²."""
        return self.a * self.a - self.a * self.b + self.b * self.b

    def __mul__(self, other: 'EisensteinInt') -> 'EisensteinInt':
        """(a + bω)(c + dω) = (ac - bd) + (ad + bc - bd)ω."""
        a, b = self.a, self.b
        c, d = other.a, other.b
        return EisensteinInt(a * c - b * d, a * d + b * c - b * d)

    def __add__(self, other: 'EisensteinInt') -> 'EisensteinInt':
        return EisensteinInt(self.a + other.a, self.b + other.b)

    def __repr__(self):
        if self.b == 0:
            return f"{self.a}"
        if self.a == 0:
            return f"{self.b}ω"
        sign = "+" if self.b > 0 else "-"
        return f"{self.a} {sign} {abs(self.b)}ω"

    def to_complex(self) -> complex:
        """Convert to complex number."""
        return self.a + self.b * complex(OMEGA_REAL, OMEGA_IMAG)

    def conjugate(self) -> 'EisensteinInt':
        """Eisenstein conjugate: a + bω → a + bω̄ = (a-b) + (-b)ω̄... 
        Actually: conj(a + bω) = a + bω² = a + b(-1-ω) = (a-b) - bω."""
        return EisensteinInt(self.a - self.b, -self.b)


# ============================================================================
# Cubic Residue Symbol
# ============================================================================

def is_prime_1_mod_3(p: int) -> bool:
    """Check if p is a prime ≡ 1 (mod 3). Only these primes split in Z[ω]."""
    if p < 2:
        return False
    if p == 3:
        return False  # 3 ramifies in Z[ω]
    if p % 3 != 1:
        return False
    # Primality check
    if p < 4:
        return p >= 2
    if p % 2 == 0:
        return False
    for i in range(3, int(math.isqrt(p)) + 1, 2):
        if p % i == 0:
            return False
    return True


def cubic_residue_symbol(a: int, p: int) -> int:
    """
    Compute the cubic residue symbol (a/p)₃ for p ≡ 1 (mod 3).

    (a/p)₃ = a^{(p-1)/3} mod p

    Returns:
        1 if a is a cubic residue mod p
        ω-power index (0, 1, 2) indicating which cube root of unity
        For computation: returns a^{(p-1)/3} mod p
    """
    if not is_prime_1_mod_3(p):
        return 0  # Undefined for p ≢ 1 mod 3

    if a % p == 0:
        return 0  # a ≡ 0

    result = pow(a, (p - 1) // 3, p)
    return result


def cubic_residue_class(a: int, p: int) -> int:
    """
    Classify the cubic residue: 0, 1, or 2.

    0 = cubic residue (a^{(p-1)/3} ≡ 1 mod p)
    1 = first non-residue (a^{(p-1)/3} ≡ g mod p where g is a primitive cube root of 1)
    2 = second non-residue (a^{(p-1)/3} ≡ g² mod p)
    """
    r = cubic_residue_symbol(a, p)
    if r == 0:
        return -1  # degenerate
    if r == 1:
        return 0  # cubic residue
    # Find which non-residue class
    # g = primitive cube root of unity mod p
    # g satisfies g³ ≡ 1, g ≠ 1
    # g = r or g = r² mod p... just check
    if pow(r, 3, p) == 1 and r != 1:
        return 1  # First cube root of unity
    r2 = pow(r, 2, p)
    if pow(r2, 3, p) == 1 and r2 != 1:
        return 2
    return 1  # Default to first non-residue


def has_cubic_asymmetry(p: int, q: int) -> bool:
    """
    Check if the cubic residue symbols (p/q)₃ and (q/p)₃ differ.

    The cubic reciprocity law (over Eisenstein integers) governs this.
    Asymmetry occurs when the cubic symbols are different cube roots of unity.
    """
    if not (is_prime_1_mod_3(p) and is_prime_1_mod_3(q)):
        return False

    c_pq = cubic_residue_symbol(p, q)
    c_qp = cubic_residue_symbol(q, p)

    return c_pq != c_qp


# ============================================================================
# 3-Adic Classification of Primes
# ============================================================================

def classify_primes_mod9(primes: List[int]) -> Dict[int, List[int]]:
    """
    Classify primes by their mod-9 residue class.

    The mod-9 structure should govern the 3-adic correction tower.
    Primes ≡ 1 mod 3 (i.e., mod 9 ∈ {1, 4, 7}) split in Z[ω].
    Primes ≡ 2 mod 3 (i.e., mod 9 ∈ {2, 5, 8}) remain inert in Z[ω].

    Returns:
        Dict mapping mod-9 residue -> list of primes
    """
    classes: Dict[int, List[int]] = {r: [] for r in range(9)}
    for p in primes:
        if p > 3:
            classes[p % 9].append(p)
    return classes


def classify_primes_mod3(primes: List[int]) -> Tuple[List[int], List[int]]:
    """
    Split primes into:
    - Split primes: p ≡ 1 (mod 3) — split in Z[ω], have cubic structure
    - Inert primes: p ≡ 2 (mod 3) — remain inert in Z[ω]

    Returns:
        (split_primes, inert_primes)
    """
    split = [p for p in primes if p > 3 and p % 3 == 1]
    inert = [p for p in primes if p > 3 and p % 3 == 2]
    return split, inert


# ============================================================================
# Combined 2-adic × 3-adic Classification
# ============================================================================

def combined_classification(primes: List[int]) -> Dict[Tuple[int, int], List[int]]:
    """
    Classify primes by BOTH mod-4 (Legendre/2-adic) and mod-3 (cubic/3-adic).

    This gives the mod-12 classification (lcm(4,3) = 12):
    - (1 mod 4, 1 mod 3) = 1 mod 12: Q₁ + Split → symmetric Legendre, has cubic structure
    - (3 mod 4, 1 mod 3) = 7 mod 12: Q₃ + Split → asymmetric Legendre, has cubic structure
    - (1 mod 4, 2 mod 3) = 5 mod 12: Q₁ + Inert → symmetric Legendre, no cubic structure
    - (3 mod 4, 2 mod 3) = 11 mod 12: Q₃ + Inert → asymmetric Legendre, no cubic structure

    Returns:
        Dict mapping (mod4_class, mod3_class) -> list of primes
    """
    classes: Dict[Tuple[int, int], List[int]] = {}
    for mod4 in [1, 3]:
        for mod3 in [1, 2]:
            classes[(mod4, mod3)] = []

    for p in primes:
        if p <= 3:
            continue
        mod4 = p % 4
        mod3 = p % 3
        if mod3 == 0:
            continue  # p = 3 special case
        key = (mod4, mod3)
        if key in classes:
            classes[key].append(p)

    return classes


def build_cubic_interaction_matrix(primes: List[int]) -> np.ndarray:
    """
    Build the cubic residue interaction matrix for primes ≡ 1 (mod 3).

    C[i,j] = cubic_residue_symbol(p_i, p_j) for primes that split in Z[ω].

    Returns:
        Complex matrix where entries are cube roots of unity (1, ω, ω²)
        mapped to {1, exp(2πi/3), exp(4πi/3)}.
    """
    split_primes = [p for p in primes if is_prime_1_mod_3(p)]
    n = len(split_primes)

    C = np.zeros((n, n), dtype=np.complex128)
    omega = np.exp(2j * np.pi / 3)

    for i in range(n):
        for j in range(n):
            if i != j:
                cls = cubic_residue_class(split_primes[i], split_primes[j])
                if cls == 0:
                    C[i, j] = 1.0
                elif cls == 1:
                    C[i, j] = omega
                elif cls == 2:
                    C[i, j] = omega**2

    return C, split_primes


def cubic_asymmetry_count(primes: List[int]) -> Dict[str, int]:
    """
    Count cubic residue asymmetries among prime pairs.

    Returns:
        Dict with 'total_pairs', 'asymmetric_pairs', 'symmetric_pairs'
    """
    split_primes = [p for p in primes if is_prime_1_mod_3(p)]
    total = 0
    asymmetric = 0
    symmetric = 0

    for i, p in enumerate(split_primes):
        for q in split_primes[i + 1:]:
            total += 1
            if has_cubic_asymmetry(p, q):
                asymmetric += 1
            else:
                symmetric += 1

    return {
        'total_pairs': total,
        'asymmetric_pairs': asymmetric,
        'symmetric_pairs': symmetric,
        'asymmetry_rate': asymmetric / total if total > 0 else 0,
    }


def correction_tower_analysis(primes: List[int]) -> Dict[str, object]:
    """
    Analyze whether the combined mod-12 classification generates the
    {0, ±3, ±9, ±27} correction tower.

    The hypothesis:
    - mod-4 structure (Legendre/quadratic) → 2-adic component → 108 = 2²×3³
    - mod-3 structure (cubic) → 3-adic component → corrections {3⁰, 3¹, 3², 3³}
    - Combined mod-12 → full mass formula m/m_e = n×108 + correction

    Returns:
        Dict with classification counts and correction pattern analysis.
    """
    combined = combined_classification(primes)

    # Count primes in each mod-12 class
    class_counts = {k: len(v) for k, v in combined.items()}

    # The correction hierarchy from particle masses:
    correction_particles = {
        0: ['proton'],       # 17 × 108 + 0
        3: ['neutron'],      # 17 × 108 + 3
        -9: ['muon'],        # 2 × 108 - 9
        9: ['charm'],        # 23 × 108 + 9
        21: ['tau'],         # 32 × 108 + 21  (21 = 3 × 7)
        -27: ['higgs', 'top', 'bottom'],  # various × 108 - 27
        27: ['W', 'Z'],      # various × 108 + 27
    }

    # Map corrections to 3-adic classification
    correction_mod3 = {}
    for c, particles in correction_particles.items():
        correction_mod3[c] = {
            'particles': particles,
            'c_mod_3': c % 3,
            'c_mod_9': c % 9,
            'is_power_of_3': c != 0 and (c in [1, 3, 9, 27] or c in [-1, -3, -9, -27]),
        }

    return {
        'class_counts': class_counts,
        'correction_mod3_analysis': correction_mod3,
        'total_primes_analyzed': sum(class_counts.values()),
    }
