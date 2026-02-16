"""
Periodicity detection in prime resonance structures.

Analyzes the 108-periodicity of resonant pair density and
the period-24 digital root cycling in Legendre symbol statistics.
"""

from typing import Dict, List, Tuple
import numpy as np
from .prime_utils import (
    sieve_primes, legendre_symbol, has_legendre_asymmetry,
    mod30_residue, digital_root, COPRIME_RESIDUES_MOD30, TWIST_UNIT
)


def cumulative_asymmetric_pair_count(primes: List[int]) -> List[Tuple[int, int]]:
    """
    Compute the cumulative count of asymmetric Legendre pairs up to each prime.

    For each prime p_n, count how many pairs (p_i, p_j) with i < j ≤ n
    have Legendre asymmetry.

    Returns:
        List of (prime, cumulative_count) pairs
    """
    working = [p for p in primes if p > 5]
    counts = []
    total = 0

    for idx, p_n in enumerate(working):
        # Count new asymmetric pairs involving p_n
        for i in range(idx):
            if has_legendre_asymmetry(working[i], p_n):
                total += 1
        counts.append((p_n, total))

    return counts


def detrended_pair_density(
    cumulative: List[Tuple[int, int]]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the detrended fluctuation of cumulative asymmetric pair counts.

    Removes the linear trend to expose periodic modulation.

    Returns:
        (prime_values, detrended_counts) as numpy arrays
    """
    if not cumulative:
        return np.array([]), np.array([])

    primes_arr = np.array([c[0] for c in cumulative], dtype=np.float64)
    counts_arr = np.array([c[1] for c in cumulative], dtype=np.float64)

    # Fit linear trend
    coeffs = np.polyfit(primes_arr, counts_arr, 1)
    trend = np.polyval(coeffs, primes_arr)

    detrended = counts_arr - trend
    return primes_arr, detrended


def detect_periodicity_fft(
    signal: np.ndarray,
    sample_spacing: float = 1.0
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Detect dominant periodicity in a signal using FFT.

    Returns:
        (frequencies, amplitudes, dominant_period)
    """
    if len(signal) < 4:
        return np.array([]), np.array([]), 0.0

    # Windowed FFT to reduce spectral leakage
    window = np.hanning(len(signal))
    windowed = signal * window

    fft_vals = np.fft.rfft(windowed)
    freqs = np.fft.rfftfreq(len(signal), d=sample_spacing)
    amplitudes = np.abs(fft_vals)

    # Exclude DC component
    if len(amplitudes) > 1:
        amplitudes[0] = 0

    if np.max(amplitudes) == 0:
        return freqs, amplitudes, 0.0

    peak_idx = np.argmax(amplitudes[1:]) + 1
    dominant_freq = freqs[peak_idx]

    dominant_period = 1.0 / dominant_freq if dominant_freq > 0 else 0.0

    return freqs, amplitudes, dominant_period


def analyze_108_periodicity(primes: List[int], tolerance: float = 0.3) -> Dict[str, object]:
    """
    Test whether the density of asymmetric Legendre pairs shows
    108-periodicity (as predicted by the twist framework where 108 = 2²×3³).

    The test:
    1. Compute cumulative asymmetric pair counts
    2. Detrend to remove linear growth
    3. FFT to find dominant period
    4. Check if dominant period is within tolerance of 108

    Returns:
        Dict with 'dominant_period', 'period_108_amplitude',
        'has_108_signal', 'frequency_spectrum'
    """
    cumulative = cumulative_asymmetric_pair_count(primes)
    primes_arr, detrended = detrended_pair_density(cumulative)

    if len(detrended) < 10:
        return {
            'dominant_period': 0.0,
            'period_108_amplitude': 0.0,
            'has_108_signal': False,
            'frequency_spectrum': (np.array([]), np.array([])),
        }

    # Average spacing between consecutive primes for FFT sampling
    avg_spacing = float(np.mean(np.diff(primes_arr))) if len(primes_arr) > 1 else 1.0

    freqs, amplitudes, dominant_period = detect_periodicity_fft(detrended, avg_spacing)

    # Check for 108-periodicity signal
    target_freq = 1.0 / TWIST_UNIT  # ≈ 0.00926
    freq_window = target_freq * tolerance

    # Find amplitude near 108-frequency
    mask = (freqs > target_freq - freq_window) & (freqs < target_freq + freq_window)
    period_108_amplitude = float(np.max(amplitudes[mask])) if np.any(mask) else 0.0

    # Compare to noise floor
    noise_floor = float(np.median(amplitudes[1:])) if len(amplitudes) > 1 else 0.0
    has_108_signal = period_108_amplitude > 2.0 * noise_floor

    return {
        'dominant_period': dominant_period,
        'period_108_amplitude': period_108_amplitude,
        'noise_floor': noise_floor,
        'has_108_signal': has_108_signal,
        'frequency_spectrum': (freqs, amplitudes),
    }


def legendre_symbol_running_sum(primes: List[int]) -> List[Tuple[int, int, float]]:
    """
    Compute the running sum S(n) = Σ_{k<n} (p_k / p_n) for consecutive primes.

    Returns:
        List of (prime, mod90_residue, running_sum) tuples
    """
    working = [p for p in primes if p > 5]
    results = []

    for n in range(1, len(working)):
        p_n = working[n]
        s = sum(legendre_symbol(working[k], p_n) for k in range(n))
        r90 = p_n % 90
        results.append((p_n, r90, float(s)))

    return results


def period_24_in_legendre_sums(
    running_sums: List[Tuple[int, int, float]]
) -> Dict[int, float]:
    """
    Group Legendre running sums by mod-90 residue class and compute means.

    If the period-24 structure influences Legendre symbol statistics,
    different mod-90 classes should have systematically different mean sums.

    Returns:
        Dict mapping mod-90 residue -> mean running sum
    """
    # Valid mod-90 residues for numbers coprime to 30
    mod90_groups: Dict[int, List[float]] = {}
    for prime_val, r90, s in running_sums:
        if r90 not in mod90_groups:
            mod90_groups[r90] = []
        mod90_groups[r90].append(s)

    means = {r: float(np.mean(vals)) for r, vals in mod90_groups.items() if vals}
    return means


def autocorrelation_of_legendre_sums(
    running_sums: List[Tuple[int, int, float]],
    max_lag: int = 48
) -> np.ndarray:
    """
    Compute autocorrelation of the Legendre running sum sequence.

    If period-24 structure is present, autocorrelation should peak at lag 24.

    Returns:
        1D array of autocorrelation values for lags 0..max_lag
    """
    sums = np.array([s for _, _, s in running_sums])
    if len(sums) < max_lag + 1:
        return np.array([])

    # Normalize
    sums = sums - np.mean(sums)
    var = np.var(sums)
    if var == 0:
        return np.zeros(max_lag + 1)

    acf = np.correlate(sums, sums, mode='full')
    mid = len(acf) // 2
    acf = acf[mid:mid + max_lag + 1] / (var * len(sums))

    return acf
