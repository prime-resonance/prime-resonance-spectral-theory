"""
Tests for the Periodic Table Generator and Successive Ionization Energies.

Validates:
  - successive_ionization_energies function correctness
  - calculate_ionization_transition returns valid configs
  - NIST successive IE data module
  - Periodic table generator output format and statistics
  - Successive IE trends (monotonically increasing)
  - Core-shell jumps in successive IEs
  - log(IE_k) data generation
"""

import math
import pytest
from src.atomic_physics.energy import (
    ionization_energy,
    calculate_ionization_transition,
    successive_ionization_energies,
)
from src.atomic_physics.orbitals import electron_configuration
from src.nist_successive_ie import (
    NIST_SUCCESSIVE_IE,
    ELEMENT_SYMBOLS,
    get_nist_successive_ie,
    get_max_ie_level,
)
from src.periodic_table_generator import (
    generate_first_ie_table,
    compute_ie_statistics,
    generate_successive_ie_table,
    compute_successive_ie_statistics,
    generate_log_ie_plot_data,
    format_first_ie_markdown,
    format_successive_ie_markdown,
    run_full_comparison,
)


# ──────────────────────────────────────────────────────────
# NIST Data Module Tests
# ──────────────────────────────────────────────────────────

class TestNISTSuccessiveData:
    """Verify the NIST successive IE reference data module."""

    def test_hydrogen_single_ie(self):
        """H has exactly 1 successive IE."""
        data = get_nist_successive_ie(1)
        assert data is not None
        assert len(data) == 1
        assert abs(data[0] - 13.598) < 0.01

    def test_helium_two_ies(self):
        """He has exactly 2 successive IEs."""
        data = get_nist_successive_ie(2)
        assert data is not None
        assert len(data) == 2
        assert data[1] > data[0], "IE_2 > IE_1 for He"

    def test_lithium_three_ies(self):
        """Li has 3 successive IEs."""
        data = get_nist_successive_ie(3)
        assert len(data) == 3
        # Core shell jump: IE_2 >> IE_1
        assert data[1] > 10 * data[0], "Li IE_2 should be >> IE_1 (core shell jump)"

    def test_sodium_ten_ies(self):
        """Na has 10 successive IEs."""
        data = get_nist_successive_ie(11)
        assert len(data) == 10
        # Core shell jump between IE_1 and IE_2
        assert data[1] > 5 * data[0], "Na IE_2 >> IE_1"

    def test_nist_data_coverage(self):
        """NIST data covers Z=1..36."""
        for Z in range(1, 37):
            data = get_nist_successive_ie(Z)
            assert data is not None, f"Missing NIST data for Z={Z}"
            assert len(data) >= 1

    def test_monotonically_increasing(self):
        """Each element's successive IEs should be monotonically increasing."""
        for Z in range(1, 37):
            data = get_nist_successive_ie(Z)
            for i in range(1, len(data)):
                assert data[i] > data[i - 1], (
                    f"Z={Z}: IE_{i+1} ({data[i]:.1f}) <= IE_{i} ({data[i-1]:.1f})"
                )

    def test_element_symbols(self):
        """Element symbols dict should have Z=1..36."""
        for Z in range(1, 37):
            assert Z in ELEMENT_SYMBOLS
        assert ELEMENT_SYMBOLS[1] == 'H'
        assert ELEMENT_SYMBOLS[26] == 'Fe'

    def test_get_max_ie_level(self):
        """get_max_ie_level returns correct count."""
        assert get_max_ie_level(1) == 1
        assert get_max_ie_level(2) == 2
        assert get_max_ie_level(11) == 10
        assert get_max_ie_level(200) == 0  # Not in database


# ──────────────────────────────────────────────────────────
# Successive Ionization Energy Function Tests
# ──────────────────────────────────────────────────────────

class TestSuccessiveIEFunction:
    """Test the successive_ionization_energies function."""

    def test_hydrogen_single_ie(self):
        """H should produce exactly 1 IE."""
        ies = successive_ionization_energies(1)
        assert len(ies) == 1
        assert ies[0] > 10, "H IE should be > 10 eV"

    def test_helium_two_ies(self):
        """He should produce exactly 2 IEs."""
        ies = successive_ionization_energies(2)
        assert len(ies) == 2
        assert ies[1] > ies[0], "He IE_2 > IE_1"

    def test_lithium_three_ies(self):
        """Li should produce 3 IEs."""
        ies = successive_ionization_energies(3, max_k=3)
        assert len(ies) == 3

    def test_max_k_limits_output(self):
        """max_k should limit the number of IEs returned."""
        ies_3 = successive_ionization_energies(10, max_k=3)
        ies_5 = successive_ionization_energies(10, max_k=5)
        assert len(ies_3) == 3
        assert len(ies_5) == 5

    def test_all_positive(self):
        """All successive IEs should be positive."""
        for Z in [1, 2, 6, 10, 11, 18, 26, 36]:
            ies = successive_ionization_energies(Z, max_k=min(Z, 10))
            for k, ie in enumerate(ies):
                assert ie > 0, f"Z={Z}, IE_{k+1} = {ie} <= 0"

    def test_monotonically_increasing(self):
        """Successive IEs should be monotonically increasing for light atoms."""
        for Z in [2, 3, 6, 7, 10, 11, 18]:
            ies = successive_ionization_energies(Z, max_k=min(Z, 8))
            for i in range(1, len(ies)):
                assert ies[i] > ies[i - 1], (
                    f"Z={Z}: IE_{i+1} ({ies[i]:.1f}) <= IE_{i} ({ies[i-1]:.1f})"
                )

    def test_first_ie_matches_ionization_energy(self):
        """IE_1 from successive should match ionization_energy()."""
        for Z in [1, 6, 11, 18, 26, 36]:
            ies = successive_ionization_energies(Z, max_k=1)
            ie_single = ionization_energy(Z)
            assert abs(ies[0] - ie_single) < 0.001, (
                f"Z={Z}: successive IE_1 ({ies[0]:.4f}) != ionization_energy ({ie_single:.4f})"
            )

    def test_core_shell_jump_sodium(self):
        """Na: IE_2 should be much larger than IE_1 (core shell jump)."""
        ies = successive_ionization_energies(11, max_k=3)
        # IE_2 / IE_1 should be at least 5x for Na
        ratio = ies[1] / ies[0]
        assert ratio > 3, f"Na IE_2/IE_1 ratio = {ratio:.1f}, expected > 3"

    def test_core_shell_jump_magnesium(self):
        """Mg: IE_3 should be much larger than IE_2 (core shell jump)."""
        ies = successive_ionization_energies(12, max_k=4)
        ratio = ies[2] / ies[1]
        assert ratio > 3, f"Mg IE_3/IE_2 ratio = {ratio:.1f}, expected > 3"

    def test_carbon_six_ies(self):
        """C should produce 6 IEs (6 electrons)."""
        ies = successive_ionization_energies(6, max_k=6)
        assert len(ies) == 6

    def test_max_k_exceeds_electrons(self):
        """max_k > Z should produce exactly Z IEs."""
        ies = successive_ionization_energies(3, max_k=100)
        assert len(ies) == 3


# ──────────────────────────────────────────────────────────
# calculate_ionization_transition Tests
# ──────────────────────────────────────────────────────────

class TestIonizationTransition:
    """Test the transition function returns correct configs."""

    def test_hydrogen_transition(self):
        """H: removing 1s1 -> empty config."""
        config = electron_configuration(1)
        ie, new_config = calculate_ionization_transition(1, config)
        assert ie > 0
        assert new_config == [] or sum(c for _, _, c in new_config) == 0

    def test_helium_transition(self):
        """He: removing from 1s2 -> 1s1."""
        config = electron_configuration(2)
        ie, new_config = calculate_ionization_transition(2, config)
        assert ie > 0
        total_after = sum(c for _, _, c in new_config)
        assert total_after == 1

    def test_sodium_transition_removes_valence(self):
        """Na: first ionization should remove from valence 3s."""
        config = electron_configuration(11)
        ie, new_config = calculate_ionization_transition(11, config)
        # Na neutral has 11 electrons, ion has 10
        total_after = sum(c for _, _, c in new_config)
        assert total_after == 10
        assert ie > 0

    def test_iron_transition(self):
        """Fe: first ionization should produce 25-electron config."""
        config = electron_configuration(26)
        ie, new_config = calculate_ionization_transition(26, config)
        total_after = sum(c for _, _, c in new_config)
        assert total_after == 25
        assert ie > 0

    def test_empty_config(self):
        """Empty config returns 0 IE."""
        ie, config = calculate_ionization_transition(1, [])
        assert ie == 0.0
        assert config == []


# ──────────────────────────────────────────────────────────
# Periodic Table Generator Tests
# ──────────────────────────────────────────────────────────

class TestFirstIETable:
    """Test first ionization energy table generation."""

    def test_table_length(self):
        """Table should have z_max entries."""
        results = generate_first_ie_table(z_max=10)
        assert len(results) == 10

    def test_table_keys(self):
        """Each row should have required keys."""
        results = generate_first_ie_table(z_max=5)
        required_keys = {'Z', 'symbol', 'config', 'predicted_IE', 'nist_IE', 'error_pct'}
        for row in results:
            assert required_keys.issubset(row.keys())

    def test_table_z_sequence(self):
        """Z values should be sequential 1..z_max."""
        results = generate_first_ie_table(z_max=20)
        for i, row in enumerate(results):
            assert row['Z'] == i + 1

    def test_positive_predictions(self):
        """All predicted IEs should be positive."""
        results = generate_first_ie_table(z_max=36)
        for row in results:
            assert row['predicted_IE'] > 0, f"Z={row['Z']} has non-positive IE"


class TestIEStatistics:
    """Test error statistics computation."""

    def test_stats_keys(self):
        """Statistics dict should have required keys."""
        results = generate_first_ie_table(z_max=10)
        stats = compute_ie_statistics(results)
        assert 'mape' in stats
        assert 'median' in stats
        assert 'max_error' in stats
        assert 'count' in stats
        assert 'within_5pct' in stats

    def test_mape_positive(self):
        """MAPE should be positive."""
        results = generate_first_ie_table(z_max=36)
        stats = compute_ie_statistics(results)
        assert stats['mape'] > 0

    def test_count_matches(self):
        """Count should match number of rows with valid NIST data."""
        results = generate_first_ie_table(z_max=36)
        stats = compute_ie_statistics(results)
        assert stats['count'] == 36


class TestSuccessiveIETable:
    """Test successive ionization energy table generation."""

    def test_table_not_empty(self):
        """Successive IE table should have entries."""
        results = generate_successive_ie_table(z_max=5, k_max=5)
        assert len(results) > 10  # 1+2+3+4+5 = 15 entries

    def test_table_keys(self):
        """Each row should have required keys."""
        results = generate_successive_ie_table(z_max=3, k_max=3)
        required_keys = {'Z', 'symbol', 'k', 'predicted_IE', 'nist_IE', 'error_pct'}
        for row in results:
            assert required_keys.issubset(row.keys())

    def test_successive_stats(self):
        """Statistics should include overall and per-level data."""
        results = generate_successive_ie_table(z_max=10, k_max=5)
        stats = compute_successive_ie_statistics(results)
        assert 'overall' in stats
        assert 'by_level' in stats
        assert 1 in stats['by_level']  # Level 1 should exist


class TestLogIEPlotData:
    """Test log(IE_k) vs k plot data generation."""

    def test_plot_data_not_empty(self):
        """Plot data should contain entries."""
        data = generate_log_ie_plot_data(z_list=[1, 2, 6, 11], k_max=5)
        assert len(data) == 4

    def test_plot_data_keys(self):
        """Plot data keys should be Z values."""
        data = generate_log_ie_plot_data(z_list=[6, 11], k_max=5)
        assert 6 in data
        assert 11 in data

    def test_plot_data_format(self):
        """Each entry should be (k, log10(IE))."""
        data = generate_log_ie_plot_data(z_list=[6], k_max=6)
        for k, log_ie in data[6]:
            assert isinstance(k, int)
            assert isinstance(log_ie, float)
            assert k >= 1


class TestMarkdownFormatting:
    """Test markdown output formatting."""

    def test_first_ie_markdown(self):
        """First IE markdown should contain header and table."""
        results = generate_first_ie_table(z_max=5)
        stats = compute_ie_statistics(results)
        md = format_first_ie_markdown(results, stats)
        assert "## First Ionization Energy" in md
        assert "MAPE" in md
        assert "| Z |" in md

    def test_successive_ie_markdown(self):
        """Successive IE markdown should contain header and table."""
        results = generate_successive_ie_table(z_max=5, k_max=3)
        stats = compute_successive_ie_statistics(results)
        md = format_successive_ie_markdown(results, stats)
        assert "## Successive Ionization Energies" in md
        assert "Per-Level MAPE" in md


class TestFullComparison:
    """Test the run_full_comparison function."""

    def test_full_comparison_runs(self):
        """Full comparison should complete without errors."""
        output = run_full_comparison(z_max_first=10, z_max_successive=5, k_max=3)
        assert 'first_ie' in output
        assert 'successive_ie' in output
        assert 'plot_data' in output

    def test_full_comparison_has_results(self):
        """Full comparison should contain non-empty results."""
        output = run_full_comparison(z_max_first=10, z_max_successive=5, k_max=3)
        assert len(output['first_ie']['results']) == 10
        assert len(output['successive_ie']['results']) > 0
        assert len(output['plot_data']) > 0

    def test_full_comparison_has_stats(self):
        """Full comparison should contain statistics."""
        output = run_full_comparison(z_max_first=10, z_max_successive=5, k_max=3)
        assert output['first_ie']['stats']['count'] > 0
        assert output['successive_ie']['stats']['overall']['count'] > 0

    def test_full_comparison_has_markdown(self):
        """Full comparison should contain formatted markdown."""
        output = run_full_comparison(z_max_first=5, z_max_successive=3, k_max=2)
        assert len(output['first_ie']['markdown']) > 100
        assert len(output['successive_ie']['markdown']) > 100


# ──────────────────────────────────────────────────────────
# NIST Comparison Accuracy Tests
# ──────────────────────────────────────────────────────────

class TestNISTComparisonAccuracy:
    """Verify that framework predictions are within reasonable bounds of NIST."""

    def test_hydrogen_ie1_accuracy(self):
        """H IE_1 should be within 10% of NIST."""
        ies = successive_ionization_energies(1)
        nist = NIST_SUCCESSIVE_IE[1][0]
        error = abs(ies[0] - nist) / nist * 100
        assert error < 10, f"H IE_1 error = {error:.1f}%"

    def test_helium_ie1_ie2_accuracy(self):
        """He IE_1 and IE_2 should be within 15% of NIST."""
        ies = successive_ionization_energies(2, max_k=2)
        for k in range(2):
            nist = NIST_SUCCESSIVE_IE[2][k]
            error = abs(ies[k] - nist) / nist * 100
            assert error < 15, f"He IE_{k+1} error = {error:.1f}%"

    def test_overall_successive_mape_z1_36(self):
        """Overall successive IE MAPE for Z=1..36 should be < 50%."""
        results = generate_successive_ie_table(z_max=36, k_max=10)
        stats = compute_successive_ie_statistics(results)
        mape = stats['overall']['mape']
        assert mape < 50, f"Overall successive MAPE = {mape:.1f}%, expected < 50%"

    def test_first_ie_mape_z1_36(self):
        """First IE MAPE for Z=1..36 should be < 20%."""
        results = generate_first_ie_table(z_max=36)
        stats = compute_ie_statistics(results)
        assert stats['mape'] < 20, f"First IE MAPE = {stats['mape']:.1f}%, expected < 20%"
