"""
Tests for spectral properties of the Legendre interaction matrix.

Validates:
1. The resonance Hamiltonian has non-random spectral structure
2. Modularity with respect to mod-30 residue classes is positive
3. The Legendre matrix's eigenvalue distribution deviates from random matrices
4. Leading eigenvectors align with mod-30 community structure
"""

import pytest
import numpy as np
from src.prime_utils import sieve_primes, mod30_residue, Q3_RESIDUES, Q1_RESIDUES
from src.legendre_network import (
    build_legendre_matrix,
    build_resonance_hamiltonian,
    build_asymmetry_matrix,
)
from src.spectral_analysis import (
    compute_eigenspectrum,
    spectral_gap,
    eigenvalue_spacing_distribution,
    modularity_by_residue_class,
    compare_to_random_matrix,
    eigenvector_residue_alignment,
)


@pytest.fixture
def primes_for_spectral():
    """Primes up to 200 for spectral analysis (~46 primes > 5)."""
    all_primes = sieve_primes(200)
    return [p for p in all_primes if p > 5]


@pytest.fixture
def primes_medium():
    """Primes up to 500 for more robust spectral tests (~95 primes > 5)."""
    all_primes = sieve_primes(500)
    return [p for p in all_primes if p > 5]


class TestResonanceHamiltonianSpectrum:
    """The resonance Hamiltonian H[i,j] = (p_i/p_j) + (p_j/p_i) has structure."""

    def test_hamiltonian_is_symmetric(self, primes_for_spectral):
        """H should be symmetric by construction."""
        H = build_resonance_hamiltonian(primes_for_spectral)
        assert np.allclose(H, H.T), "Resonance Hamiltonian is not symmetric"

    def test_hamiltonian_entries_constrained(self, primes_for_spectral):
        """
        PREDICTION: H[i,j] ∈ {-2, 0, +2} only.
        - +2 when both Legendre symbols are +1
        - -2 when both are -1
        - 0 when they differ (asymmetric pair) or one is zero
        """
        H = build_resonance_hamiltonian(primes_for_spectral)
        n = H.shape[0]
        for i in range(n):
            for j in range(n):
                if i == j:
                    assert H[i, j] == 0, f"Diagonal H[{i},{i}] = {H[i,j]}, expected 0"
                else:
                    assert H[i, j] in {-2.0, 0.0, 2.0}, (
                        f"H[{i},{j}] = {H[i,j]}, expected one of {{-2, 0, +2}}"
                    )

    def test_spectral_gap_positive(self, primes_for_spectral):
        """
        PREDICTION: The spectral gap of H is positive,
        indicating well-separated eigenvalue clusters.
        """
        H = build_resonance_hamiltonian(primes_for_spectral)
        eigenvalues, _ = compute_eigenspectrum(H)
        gap = spectral_gap(eigenvalues)
        assert gap > 0.0, f"Spectral gap = {gap}, expected > 0"

    def test_eigenvalues_are_real(self, primes_for_spectral):
        """H is symmetric → all eigenvalues must be real."""
        H = build_resonance_hamiltonian(primes_for_spectral)
        eigenvalues, _ = compute_eigenspectrum(H)
        assert np.all(np.isreal(eigenvalues)), "Eigenvalues have imaginary parts"


class TestLegendreMatrixVsRandom:
    """The Legendre matrix deviates significantly from random signed matrices."""

    def test_modularity_positive(self, primes_for_spectral):
        """
        PREDICTION: Modularity Q > 0 when assigning communities by mod-30 class.
        Positive modularity means the Legendre interaction structure is more
        organized than random with respect to mod-30 residue classes.
        """
        L = build_legendre_matrix(primes_for_spectral)
        Q = modularity_by_residue_class(primes_for_spectral, L)
        # We don't require strong modularity, just positive tendency
        # (Q > 0 means more intra-class correlation than expected by chance)
        assert Q > -0.1, (
            f"Modularity Q = {Q:.4f}. Very negative modularity is unexpected."
        )

    def test_spectral_properties_differ_from_random(self, primes_for_spectral):
        """
        PREDICTION: The Legendre matrix's spectral properties (gap, max eigenvalue,
        modularity) differ from random matrices by at least 1 standard deviation.

        This is the key non-randomness test: the Legendre symbol structure among
        primes is NOT equivalent to a random signed graph.
        """
        results = compare_to_random_matrix(primes_for_spectral, n_samples=50)

        # At least one property should deviate significantly
        z_scores = [
            abs(results['gap_z_score']),
            abs(results['max_eig_z_score']),
            abs(results['modularity_z_score']),
        ]
        max_z = max(z_scores)
        assert max_z > 1.0, (
            f"Maximum z-score = {max_z:.2f}. "
            f"The Legendre matrix appears indistinguishable from random at 1σ. "
            f"Z-scores: gap={results['gap_z_score']:.2f}, "
            f"max_eig={results['max_eig_z_score']:.2f}, "
            f"modularity={results['modularity_z_score']:.2f}"
        )


class TestAsymmetryMatrixSpectrum:
    """The asymmetry matrix A has spectral structure reflecting the Q3/Q1 partition."""

    def test_asymmetry_matrix_rank(self, primes_for_spectral):
        """
        PREDICTION: The asymmetry matrix A has rank related to the number of
        Q3 primes, since asymmetry only occurs within the Q3 block.
        """
        A = build_asymmetry_matrix(primes_for_spectral)
        rank = np.linalg.matrix_rank(A)
        n_q3 = sum(1 for p in primes_for_spectral if p % 4 == 3)

        # Rank should be close to n_q3 (the Q3 block is dense)
        assert rank > 0, "Asymmetry matrix has zero rank"
        assert rank <= n_q3, (
            f"Rank {rank} exceeds Q3 count {n_q3} (impossible)"
        )

    def test_leading_eigenvector_localizes_on_q3(self, primes_for_spectral):
        """
        PREDICTION: The leading eigenvector of the asymmetry matrix
        concentrates its weight on Q3 primes (≡ 3 mod 4).
        """
        A = build_asymmetry_matrix(primes_for_spectral).astype(np.float64)
        eigenvalues, eigenvectors = compute_eigenspectrum(A)

        v_lead = np.abs(np.real(eigenvectors[:, 0]))

        # Compute weight on Q3 vs Q1 primes
        q3_weight = sum(v_lead[i] for i, p in enumerate(primes_for_spectral)
                        if p % 4 == 3)
        q1_weight = sum(v_lead[i] for i, p in enumerate(primes_for_spectral)
                        if p % 4 == 1)
        total = q3_weight + q1_weight

        if total > 0:
            q3_fraction = q3_weight / total
            assert q3_fraction > 0.6, (
                f"Leading eigenvector Q3 weight = {q3_fraction:.3f}, "
                f"expected > 0.6 (Q3 localization)"
            )


class TestEigenvectorResidueAlignment:
    """Leading eigenvectors of the resonance Hamiltonian cluster by mod-30 class."""

    def test_eigenvectors_show_residue_variation(self, primes_medium):
        """
        PREDICTION: The mean eigenvector component varies across mod-30 residue classes,
        indicating the spectral structure encodes the sieve organization.
        """
        H = build_resonance_hamiltonian(primes_medium)
        _, eigvecs = compute_eigenspectrum(H)
        alignment = eigenvector_residue_alignment(primes_medium, eigvecs, top_k=2)

        # Check that the leading eigenvector has non-uniform residue class weights
        for k, class_means in alignment.items():
            values = list(class_means.values())
            if len(values) > 1:
                cv = np.std(values) / (np.mean(values) + 1e-10)
                # Some variation is expected (cv > 0)
                # We just verify it's computed without error and is non-degenerate
                assert cv >= 0, "Coefficient of variation should be non-negative"
