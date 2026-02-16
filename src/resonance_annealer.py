"""
Resonance-based annealer for constraint satisfaction problems.

Implements the computational annealing approach from formalism_v4_resonance.md,
using topological prime qubits and Metropolis-Hastings dynamics.
"""

from typing import List, Tuple, Dict, Optional, Callable
import numpy as np
from .prime_utils import legendre_symbol, has_legendre_asymmetry, mod30_residue, Q3_RESIDUES


class TopologicalQubit:
    """
    A topological qubit represented by a prime pair (p, q).

    State encoding:
    - |0⟩: Configuration (p, q) with chirality χ = legendre_symbol(p, q)
    - |1⟩: Configuration (q, p) with chirality χ = legendre_symbol(q, p)

    For valid qubits (both p,q ≡ 3 mod 4), these chiralities differ.
    """

    def __init__(self, p: int, q: int):
        self.p = p
        self.q = q
        self.state = 0  # 0 or 1
        self._chi_0 = legendre_symbol(p, q)
        self._chi_1 = legendre_symbol(q, p)
        self.is_valid_qubit = (self._chi_0 != self._chi_1 and
                               self._chi_0 != 0 and self._chi_1 != 0)

    @property
    def chirality(self) -> int:
        """Current chirality value: +1 or -1."""
        return self._chi_0 if self.state == 0 else self._chi_1

    def flip(self):
        """Apply Pauli-X (swap p and q)."""
        self.state = 1 - self.state

    def reset(self, state: int = 0):
        """Set qubit to a specific state."""
        self.state = state

    def __repr__(self):
        return f"TQubit({self.p},{self.q}|state={self.state},χ={self.chirality})"


class ResonanceCSPSolver:
    """
    Solve constraint satisfaction problems using prime resonance annealing.

    The solver maps CSP constraints to a resonance Hamiltonian over
    topological qubits and uses Metropolis-Hastings annealing to find
    the ground state (maximum resonance = all constraints satisfied).
    """

    def __init__(
        self,
        qubits: List[TopologicalQubit],
        constraints: List[Tuple[int, int, int]],  # (qubit_i, qubit_j, coupling)
        biases: Optional[Dict[int, int]] = None,   # qubit_index -> preferred chirality
    ):
        """
        Args:
            qubits: List of TopologicalQubit instances
            constraints: List of (i, j, J) where J=+1 means "match chirality",
                        J=-1 means "opposite chirality"
            biases: Optional dict of single-qubit biases (h_i values)
        """
        self.qubits = qubits
        self.constraints = constraints
        self.biases = biases or {}
        self.n_qubits = len(qubits)

    def resonance_energy(self) -> float:
        """
        Compute the total resonance energy (Hamiltonian).

        R = Σ_k J_k * χ_i * χ_j + Σ_i h_i * χ_i

        Higher resonance = more constraints satisfied.
        """
        energy = 0.0

        # Pairwise interactions
        for i, j, J in self.constraints:
            chi_i = self.qubits[i].chirality
            chi_j = self.qubits[j].chirality
            energy += J * chi_i * chi_j

        # Single-qubit biases
        for i, h in self.biases.items():
            energy += h * self.qubits[i].chirality

        return energy

    def max_possible_energy(self) -> float:
        """Compute the maximum possible energy (all constraints satisfied)."""
        return float(sum(abs(J) for _, _, J in self.constraints) +
                     sum(abs(h) for h in self.biases.values()))

    def anneal(
        self,
        n_steps: int = 1000,
        T_start: float = 5.0,
        T_end: float = 0.01,
        cooling: str = "exponential",
        seed: Optional[int] = None,
    ) -> Dict[str, object]:
        """
        Run Metropolis-Hastings annealing.

        Args:
            n_steps: Number of annealing steps
            T_start: Initial temperature
            T_end: Final temperature
            cooling: Cooling schedule ("exponential" or "linear")
            seed: Random seed for reproducibility

        Returns:
            Dict with 'final_energy', 'max_energy', 'trajectory', 'converged',
            'convergence_step', 'final_states'
        """
        rng = np.random.default_rng(seed)

        # Initialize random states
        for q in self.qubits:
            q.reset(rng.integers(0, 2))

        # Cooling schedule
        if cooling == "exponential":
            decay = (T_end / T_start) ** (1.0 / max(n_steps - 1, 1))
            temperatures = [T_start * (decay ** step) for step in range(n_steps)]
        else:  # linear
            temperatures = [T_start + (T_end - T_start) * step / max(n_steps - 1, 1)
                          for step in range(n_steps)]

        trajectory = []
        max_energy = self.max_possible_energy()
        best_energy = -float('inf')
        convergence_step = -1

        for step, T in enumerate(temperatures):
            # Propose a random single-qubit flip
            flip_idx = rng.integers(0, self.n_qubits)

            # Current energy
            E_current = self.resonance_energy()

            # Flip and compute new energy
            self.qubits[flip_idx].flip()
            E_proposed = self.resonance_energy()

            # Metropolis criterion (we MAXIMIZE resonance, so accept if ΔE > 0)
            delta_E = E_proposed - E_current
            if delta_E > 0:
                accept = True
            else:
                accept = rng.random() < np.exp(delta_E / max(T, 1e-10))

            if not accept:
                self.qubits[flip_idx].flip()  # Undo
                E_final = E_current
            else:
                E_final = E_proposed

            trajectory.append(E_final)

            if E_final > best_energy:
                best_energy = E_final

            if E_final >= max_energy and convergence_step < 0:
                convergence_step = step

        final_energy = self.resonance_energy()
        converged = (final_energy >= max_energy)

        return {
            'final_energy': final_energy,
            'max_energy': max_energy,
            'trajectory': trajectory,
            'converged': converged,
            'convergence_step': convergence_step,
            'final_states': [q.state for q in self.qubits],
        }


def build_2sat_problem(
    prime_pairs: List[Tuple[int, int]],
    constraints: List[Tuple[int, int, int]],
    biases: Optional[Dict[int, int]] = None,
) -> ResonanceCSPSolver:
    """
    Build a 2-SAT problem instance using the given prime pairs as qubits.

    Args:
        prime_pairs: List of (p, q) prime pairs for qubits
        constraints: List of (i, j, J) coupling constraints
        biases: Optional single-qubit biases

    Returns:
        ResonanceCSPSolver instance
    """
    qubits = [TopologicalQubit(p, q) for p, q in prime_pairs]
    return ResonanceCSPSolver(qubits, constraints, biases)


def select_qubit_primes_mod30_guided(
    primes: List[int],
    n_qubits: int,
    seed: Optional[int] = None
) -> List[Tuple[int, int]]:
    """
    Select prime pairs for qubits using mod-30 guidance (Q3 class only).

    Only selects primes ≡ 3 (mod 4), which are guaranteed to produce
    Legendre asymmetry. Groups them by mod-30 residue class for
    intra-class coherence.

    Returns:
        List of (p, q) pairs suitable for topological qubits
    """
    rng = np.random.default_rng(seed)

    # Filter to Q3 primes (≡ 3 mod 4)
    q3_primes = [p for p in primes if p > 5 and p % 4 == 3]

    if len(q3_primes) < 2 * n_qubits:
        raise ValueError(f"Not enough Q3 primes. Have {len(q3_primes)}, need {2 * n_qubits}")

    # Group by mod-30 residue class
    groups: Dict[int, List[int]] = {r: [] for r in Q3_RESIDUES}
    for p in q3_primes:
        r = mod30_residue(p)
        if r in Q3_RESIDUES:
            groups[r].append(p)

    pairs = []
    available = list(q3_primes)
    rng.shuffle(available)

    for i in range(0, min(2 * n_qubits, len(available)), 2):
        if i + 1 < len(available):
            pairs.append((available[i], available[i + 1]))

    return pairs[:n_qubits]


def select_qubit_primes_random(
    primes: List[int],
    n_qubits: int,
    seed: Optional[int] = None
) -> List[Tuple[int, int]]:
    """
    Select prime pairs randomly (no mod-30 guidance).
    Some pairs may lack Legendre asymmetry.

    Returns:
        List of (p, q) pairs (not all may be valid qubits)
    """
    rng = np.random.default_rng(seed)
    working = [p for p in primes if p > 5]

    if len(working) < 2 * n_qubits:
        raise ValueError(f"Not enough primes. Have {len(working)}, need {2 * n_qubits}")

    rng.shuffle(working)
    pairs = []
    for i in range(0, 2 * n_qubits, 2):
        pairs.append((working[i], working[i + 1]))

    return pairs[:n_qubits]
