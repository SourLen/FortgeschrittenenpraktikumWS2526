from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


# -------------------------
# I/O helpers
# -------------------------

def load_tka_counts(path: Path, skip_header_lines: int = 2) -> np.ndarray:
    """
    Load .tka file as counts-per-channel array.
    Assumes first `skip_header_lines` lines are header/non-count data.
    """
    data = np.loadtxt(path)
    return np.asarray(data[skip_header_lines:], dtype=float)


# -------------------------
# Data model
# -------------------------

@dataclass(frozen=True)
class Spectrum:
    name: str
    counts: np.ndarray
    live_time_s: float  # measurement live time in seconds

    @property
    def channel(self) -> np.ndarray:
        return np.arange(len(self.counts))


@dataclass(frozen=True)
class CorrectedSpectrum:
    """
    Spectrum after subtracting a time-scaled background:
        corr = signal - alpha * background, alpha = T_signal / T_bg
    """
    name: str
    counts_corr: np.ndarray
    sigma_corr: np.ndarray
    alpha: float
    live_time_s: float

    @property
    def channel(self) -> np.ndarray:
        return np.arange(len(self.counts_corr))


def subtract_background(signal: Spectrum, background: Spectrum) -> CorrectedSpectrum:
    alpha = signal.live_time_s / background.live_time_s
    counts_corr = signal.counts - alpha * background.counts

    # Poisson error propagation for difference of independent Poisson counts
    sigma = np.sqrt(np.maximum(signal.counts, 0.0) + (alpha**2) * np.maximum(background.counts, 0.0))
    sigma = np.where(sigma == 0, 1.0, sigma)  # avoid zeros that can break later steps

    return CorrectedSpectrum(
        name=f"{signal.name} (bg-subtracted)",
        counts_corr=counts_corr,
        sigma_corr=sigma,
        alpha=alpha,
        live_time_s=signal.live_time_s,
    )


# -------------------------
# Plotting
# -------------------------

def plot_spectrum_and_noise(spectrum: Spectrum, noise: Spectrum) -> None:
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    ax1.errorbar(spectrum.channel, spectrum.counts, yerr=np.sqrt(spectrum.counts), fmt=".", label=spectrum.name)
    ax1.set_title(f"Spectrum: {spectrum.name}")
    ax1.set_ylabel("Counts")
    ax1.legend()
    ax2.errorbar(noise.channel, noise.counts, yerr=np.sqrt(noise.counts), fmt=".", label=noise.name, color="orange")
    ax2.set_title(f"Noise: {noise.name}")
    ax2.set_ylabel("Counts")
    ax2.legend()
    # Subtracted spectrum
    corr_spec = subtract_background(spectrum, noise)
    ax3.errorbar(corr_spec.channel, corr_spec.counts_corr, yerr=corr_spec.sigma_corr, fmt=".", label=corr_spec.name, color="green")
    ax3.set_title(f"Background-subtracted Spectrum: {corr_spec.name}")
    ax3.set_xlabel("Channel")
    ax3.set_ylabel("Counts (bg-subtracted)")
    ax3.legend()
    plt.tight_layout()
    plt.show()

