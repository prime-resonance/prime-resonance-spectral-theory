"""
Core prime number utilities: generation, Legendre symbols, and classification.
"""

from typing import List, Tuple
import math


def sieve_primes(n: int) -> List[int]:
    """Generate all primes up to n using Sieve of Eratosthenes."""
    if n < 2:
        return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
    return [i for i in range(2, n + 1) if is_prime[i]]


def legendre_symbol(a: int, p: int) -> int:
    """
    Compute the Legendre symbol (a/p) for odd prime p.
    For p=2, we use the Kronecker extension: (a/2) based on a mod 8.

    Returns:
        +1 if a is a quadratic residue mod p (and a ≢ 0 mod p)
        -1 if a is a quadratic non-residue mod p
         0 if a ≡ 0 mod p
    """
    if p < 2:
        raise ValueError(f"p must be a prime, got {p}")

    if p == 2:
        # Kronecker extension for p=2
        a_mod = a % 2
        if a_mod == 0:
            return 0
        a_mod8 = a % 8
        if a_mod8 in (1, 7):
            return 1
        return -1  # a_mod8 in (3, 5)

    a = a % p
    if a == 0:
        return 0

    result = pow(a, (p - 1) // 2, p)
    # Python's pow returns p-1 for -1 in modular arithmetic
    if result == p - 1:
        return -1
    return result  # Will be 0 or 1


def has_legendre_asymmetry(p: int, q: int) -> bool:
    """
    Check if the prime pair (p, q) has Legendre symbol asymmetry.

    By quadratic reciprocity:
        (p/q)(q/p) = (-1)^((p-1)/2 * (q-1)/2)

    Asymmetry occurs when (p/q) ≠ (q/p), which requires:
        (-1)^((p-1)/2 * (q-1)/2) = -1

    This happens iff both p ≡ 3 (mod 4) AND q ≡ 3 (mod 4).
    """
    ls_pq = legendre_symbol(p, q)
    ls_qp = legendre_symbol(q, p)
    return ls_pq != ls_qp and ls_pq != 0 and ls_qp != 0


def mod4_class(p: int) -> int:
    """Return the mod-4 residue class of p (1 or 3 for odd primes)."""
    return p % 4


def mod30_residue(n: int) -> int:
    """Return the mod-30 residue of n."""
    return n % 30


def digital_root(n: int) -> int:
    """
    Compute the digital root of n.
    The digital root is the single digit obtained by repeatedly summing digits.
    """
    if n == 0:
        return 0
    return 1 + (n - 1) % 9


def is_coprime_to_30(n: int) -> bool:
    """Check if n is coprime to 30 (not divisible by 2, 3, or 5)."""
    return n % 2 != 0 and n % 3 != 0 and n % 5 != 0


# The 8 coprime residues mod 30
COPRIME_RESIDUES_MOD30 = frozenset({1, 7, 11, 13, 17, 19, 23, 29})

# Partition by mod-4 class (determines quadratic reciprocity asymmetry)
Q3_RESIDUES = frozenset({7, 11, 19, 23})   # ≡ 3 (mod 4)
Q1_RESIDUES = frozenset({1, 13, 17, 29})   # ≡ 1 (mod 4)

# Period-24 digital root cycle for numbers coprime to 30
PERIOD_24_DIGITAL_ROOTS = [1, 7, 2, 4, 8, 1, 5, 2, 4, 1, 5, 7, 2, 4, 8, 5, 7, 4, 8, 1, 5, 7, 2, 8]

# The fundamental twist unit
TWIST_UNIT = 108  # = 2^2 * 3^3

# ── Full Primorial Framework (2310 = 2 × 3 × 5 × 7 × 11) ──

# Primorial hierarchy
PRIMORIAL_2 = 2
PRIMORIAL_6 = 6        # = 2 × 3
PRIMORIAL_30 = 30      # = 2 × 3 × 5
PRIMORIAL_210 = 210    # = 2 × 3 × 5 × 7
PRIMORIAL_2310 = 2310  # = 2 × 3 × 5 × 7 × 11  (full primorial)

# Twist unit hierarchy: T(P_k) = 3³ × Π_{p≥5, p|P_k} (p−1)
TWIST_P1 = 27       # T(2) = 27
TWIST_P2 = 27       # T(6) = 27
TWIST_P3 = 108      # T(30) = 27 × 4 = 108
TWIST_P4 = 648      # T(210) = 27 × 4 × 6 = 648
TWIST_P5 = 6480     # T(2310) = 27 × 4 × 6 × 10 = 6480  (= 2⁴ × 3⁴ × 5)

# Key structural constants
TREFOIL_COMPLEXITY = 17          # c(T_{3,1}) = 17
BOUNDARY_PRIME = 29              # Largest coprime residue mod 30
CORRECTION_BASE = 27             # = 3³, 3-adic correction base
PHI_P3 = 8                       # φ(30) = 8
PHI_P4 = 48                      # φ(210) = 48
PHI_P5 = 480                     # φ(2310) = 480

# Fine structure from sub-primorial 30
ALPHA_INV_BASE = 108 + 29        # = 137  (base from sub-primorial 30)
ALPHA_INV_REFINED = 108 + 29 + 1/27  # ≈ 137.037  (refined by 3-adic correction)

# Coprime residues mod 2310
def coprime_residues_mod_2310() -> frozenset:
    """Return all residues coprime to 2310 in [1, 2310]."""
    return frozenset(
        r for r in range(1, 2311)
        if r % 2 != 0 and r % 3 != 0 and r % 5 != 0
        and r % 7 != 0 and r % 11 != 0
    )

def is_coprime_to_2310(n: int) -> bool:
    """Check if n is coprime to 2310 (not divisible by 2, 3, 5, 7, or 11)."""
    return (n % 2 != 0 and n % 3 != 0 and n % 5 != 0
            and n % 7 != 0 and n % 11 != 0)

def mod2310_residue(n: int) -> int:
    """Return the mod-2310 residue of n."""
    return n % 2310
