"""
Analytic Proof of the 108 Identity and Primorial Generalization.

THEOREM: The sum of one period-24 digital root cycle of the coprime-to-30
sequence equals exactly 108 = 2² × 3³.

PROOF:
1. The 24 residues mod 90 coprime to 30 distribute across 6 digit classes
   (digital roots 1,2,4,5,7,8 — excluding 3,6,9 since coprime to 3).
2. By counting, exactly φ(90)/6 = 24/6 = 4 residues per digit class.
3. Sum of non-3-divisible single digits: 1+2+4+5+7+8 = 27 = 3³.
4. Total: 4 × 27 = 108 = 2² × 3³.

GENERALIZATION: For primorial P_k, the digital root cycle sum of the
coprime-to-P_k sequence follows a predictable formula based on the
Euler totient and digit class structure.
"""

from typing import Dict, List, Tuple
import math


def euler_totient(n: int) -> int:
    """Compute Euler's totient function φ(n)."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def digital_root(n: int) -> int:
    """Compute digital root of n."""
    if n == 0:
        return 0
    return 1 + (n - 1) % 9


def coprime_residues_mod(m: int) -> List[int]:
    """Return all residues mod m coprime to m, in [1, m]."""
    return [r for r in range(1, m + 1) if math.gcd(r, m) == 1]


def digital_root_class_distribution(m: int) -> Dict[int, List[int]]:
    """
    Distribute coprime residues mod m into digital root classes.

    Returns:
        Dict mapping digital_root -> list of residues with that digital root
    """
    residues = coprime_residues_mod(m)
    classes: Dict[int, List[int]] = {}
    for r in residues:
        dr = digital_root(r)
        if dr not in classes:
            classes[dr] = []
        classes[dr].append(r)
    return classes


def period_24_proof_components() -> Dict[str, object]:
    """
    Produce all components of the analytic proof that Σ C₂₄ = 108.

    The proof structure:
    1. Compute φ(90) = 24 (number of coprime residues mod 90)
    2. Group by digital root → 6 classes of 4 each
    3. Sum of digit class labels: 1+2+4+5+7+8 = 27 = 3³
    4. Total: 4 × 27 = 108 = 2² × 3³

    Returns:
        Dict with all intermediate results for verification.
    """
    # Step 1: The period is 90 = 3 × 30, and φ(90) = 24
    modulus = 90  # = 3 × P₃ where P₃ = 30
    phi_90 = euler_totient(modulus)
    assert phi_90 == 24

    # Step 2: Coprime residues mod 90
    residues = coprime_residues_mod(modulus)
    assert len(residues) == 24

    # Step 3: Group by digital root
    dr_classes = digital_root_class_distribution(modulus)

    # Step 4: Verify exactly 6 classes (excluding dr ∈ {3, 6, 9})
    valid_dr = {1, 2, 4, 5, 7, 8}
    assert set(dr_classes.keys()) == valid_dr

    # Step 5: Each class has exactly 4 elements
    class_sizes = {dr: len(members) for dr, members in dr_classes.items()}
    for dr, size in class_sizes.items():
        assert size == 4, f"Digital root {dr} has {size} residues, expected 4"

    # Step 6: Sum of valid digital roots = 27 = 3³
    dr_label_sum = sum(valid_dr)
    assert dr_label_sum == 27
    assert dr_label_sum == 3**3

    # Step 7: Elements per class × sum of class labels = 4 × 27 = 108
    total = 4 * 27
    assert total == 108
    assert total == 2**2 * 3**3

    # Step 8: Also verify by direct summation
    direct_sum = sum(digital_root(r) for r in residues)
    assert direct_sum == 108

    return {
        'modulus': modulus,
        'phi_modulus': phi_90,
        'residues': residues,
        'dr_classes': dr_classes,
        'class_sizes': class_sizes,
        'valid_digital_roots': sorted(valid_dr),
        'dr_label_sum': dr_label_sum,
        'elements_per_class': 4,
        'total_dr_sum': total,
        'direct_sum': direct_sum,
        'factorization': '2² × 3³',
    }


def why_27() -> Dict[str, object]:
    """
    Prove that 1+2+4+5+7+8 = 27 = 3³ is necessary.

    This is the sum of single digits NOT divisible by 3.
    Sum(1..9) = 45. Sum(3+6+9) = 18. Remainder = 45 - 18 = 27 = 3³.

    The 3-divisible digits are excluded because coprime to 30 implies
    coprime to 3, so no residue can have digital root 3, 6, or 9.
    """
    all_digits = list(range(1, 10))
    digit_sum = sum(all_digits)
    assert digit_sum == 45

    three_divisible = [d for d in all_digits if d % 3 == 0]
    assert three_divisible == [3, 6, 9]
    three_div_sum = sum(three_divisible)
    assert three_div_sum == 18

    non_three_div = [d for d in all_digits if d % 3 != 0]
    assert non_three_div == [1, 2, 4, 5, 7, 8]
    non_three_sum = sum(non_three_div)
    assert non_three_sum == 27
    assert non_three_sum == digit_sum - three_div_sum
    assert non_three_sum == 3**3

    return {
        'all_digit_sum': digit_sum,
        'three_divisible_sum': three_div_sum,
        'non_three_sum': non_three_sum,
        'non_three_sum_is_3_cubed': non_three_sum == 3**3,
        'why': "Coprime to 30 ⊃ coprime to 3 → excludes dr ∈ {3,6,9}",
    }


def why_4_per_class() -> Dict[str, object]:
    """
    Prove that each digital root class has exactly 4 representatives mod 90.

    φ(90) = 24 coprime residues distributed among 6 valid digit classes.
    By the structure of the coprime-to-30 residues and mod-9 arithmetic:
    - 90 = 9 × 10 and lcm(30, 9) = 90
    - Each mod-9 class coprime to 3 has exactly φ(90)/6 = 4 representatives
    """
    modulus = 90
    residues = coprime_residues_mod(modulus)
    dr_classes = digital_root_class_distribution(modulus)

    # The uniformity follows from:
    # 1. The Chinese Remainder Theorem: Z/90Z ≅ Z/9Z × Z/10Z
    # 2. Coprime to 30 means coprime to 2, 3, and 5
    # 3. Coprime to 3 eliminates dr classes {3, 6, 9} (3 out of 9)
    # 4. The 6 remaining dr classes get equal share of φ(90) = 24

    # Verify CRT decomposition
    # φ(9) = 6 residues coprime to 9: {1,2,4,5,7,8}
    # φ(10) = 4 residues coprime to 10: {1,3,7,9}
    # φ(90) = φ(9) × φ(10) = 6 × 4 = 24 ✓
    phi_9 = euler_totient(9)
    phi_10 = euler_totient(10)
    assert phi_9 == 6
    assert phi_10 == 4
    assert phi_9 * phi_10 == 24

    # The 6 residues mod 9 coprime to 9 ARE the 6 valid digital roots
    coprime_to_9 = coprime_residues_mod(9)
    assert coprime_to_9 == [1, 2, 4, 5, 7, 8]

    # By CRT, each mod-9 class has φ(10) = 4 representatives mod 90
    # that are also coprime to 10 (hence coprime to 2 and 5)
    return {
        'phi_9': phi_9,
        'phi_10': phi_10,
        'product': phi_9 * phi_10,
        'coprime_to_9': coprime_to_9,
        'per_class_count': phi_10,
        'crt_explanation': "Z/90Z ≅ Z/9Z × Z/10Z; coprime to 3 selects 6 of 9 mod-9 classes; each has φ(10)=4 reps coprime to 10",
    }


def generalized_twist_unit(primorial: int) -> Dict[str, object]:
    """
    Generalize the 108 identity to other primorials.

    For primorial P_k, the "twist unit" is the sum of digital roots
    of one period of the coprime-to-P_k sequence.

    The period modulus is lcm(P_k, 9) for digital root cycling.

    Returns:
        Dict with the generalized twist unit and its factorization.
    """
    # Period modulus = lcm(P_k, 9)
    period_mod = math.lcm(primorial, 9)

    # φ(period_mod) = number of coprime residues
    phi = euler_totient(period_mod)

    # Coprime residues
    residues = coprime_residues_mod(period_mod)
    assert len(residues) == phi

    # Digital root class distribution
    dr_classes = digital_root_class_distribution(period_mod)

    # Sum of digital roots = generalized twist unit
    twist_unit = sum(digital_root(r) for r in residues)

    # Class analysis
    class_sizes = {dr: len(members) for dr, members in dr_classes.items()}
    valid_drs = sorted(dr_classes.keys())
    n_classes = len(valid_drs)
    dr_label_sum = sum(valid_drs)

    # Check uniformity
    sizes = list(class_sizes.values())
    is_uniform = len(set(sizes)) == 1
    elements_per_class = sizes[0] if is_uniform else None

    return {
        'primorial': primorial,
        'period_modulus': period_mod,
        'phi': phi,
        'twist_unit': twist_unit,
        'n_classes': n_classes,
        'valid_drs': valid_drs,
        'dr_label_sum': dr_label_sum,
        'class_sizes': class_sizes,
        'is_uniform': is_uniform,
        'elements_per_class': elements_per_class,
        'formula': f"{elements_per_class} × {dr_label_sum} = {twist_unit}" if is_uniform else "non-uniform",
    }


def general_twist_formula(primorial: int, prime_factors: List[int]) -> int:
    """
    THE GENERAL TWIST FORMULA (discovered in this work):

        T(P_k) = 3³ × Π_{p | P_k, p ≥ 5} (p - 1)

    This is proven by the chain:
    1. The twist unit = (φ(lcm(P_k, 9)) / 6) × 27
    2. φ(lcm(P_k, 9)) = φ(9) × Π_{p | P_k, p ∤ 9} φ(p)
    3. For p ≥ 5, p ∤ 9, so φ(p) = p-1
    4. φ(9)/6 = 6/6 = 1 (absorbed into the base)
    5. Result: T = 27 × Π_{p ≥ 5} (p-1) = 3³ × Π_{p ≥ 5} (p-1)

    The consecutive ratios T(P_{k+1})/T(P_k) = φ(p_{k+1}) = p_{k+1} - 1.

    Args:
        primorial: The primorial value P_k
        prime_factors: The prime factors of P_k that are ≥ 5

    Returns:
        The twist unit T(P_k)
    """
    product = 1
    for p in prime_factors:
        if p >= 5:
            product *= (p - 1)
    return 27 * product


def compute_primorial_twist_series() -> List[Dict[str, object]]:
    """
    Compute the generalized twist unit for the primorial series:
    P₁=2, P₂=6, P₃=30, P₄=210, P₅=2310

    This generates a sequence of "twist units" that may predict
    physical constants at each primorial level.

    Returns:
        List of twist unit analyses for each primorial.
    """
    primorials = [2, 6, 30, 210, 2310]
    return [generalized_twist_unit(p) for p in primorials]
