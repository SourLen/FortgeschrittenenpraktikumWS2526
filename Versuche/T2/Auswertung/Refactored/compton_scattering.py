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

def plot_comparison(
    noise_10: Spectrum,
    noise_20: Spectrum,
    sig_10: Spectrum,
    sig_20: Spectrum,
    corr_10: CorrectedSpectrum,
    corr_20: CorrectedSpectrum,
    title: str = "Compton scattering spectra (Cs-137 ring)",
    show_errorbars_corrected: bool = False,
) -> None:
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(11, 14), sharex=True)

    # Noise
    ax1.plot(noise_10.channel, noise_10.counts, marker=".", linestyle="none", label=f"{noise_10.name}")
    ax1.plot(noise_20.channel, noise_20.counts, marker=".", linestyle="none", label=f"{noise_20.name}")
    ax1.set_ylabel("counts")
    ax1.legend()

    # Raw signals
    ax2.plot(sig_10.channel, sig_10.counts, marker=".", linestyle="none", label=f"{sig_10.name}")
    ax2.plot(sig_20.channel, sig_20.counts, marker=".", linestyle="none", label=f"{sig_20.name}")
    ax2.set_ylabel("counts")
    ax2.legend()

    # Corrected
    if show_errorbars_corrected:
        ax3.errorbar(
            corr_10.channel, corr_10.counts_corr, yerr=corr_10.sigma_corr,
            fmt=".", label=f"{corr_10.name} (α={corr_10.alpha:.3g})"
        )
        ax3.errorbar(
            corr_20.channel, corr_20.counts_corr, yerr=corr_20.sigma_corr,
            fmt=".", label=f"{corr_20.name} (α={corr_20.alpha:.3g})"
        )
    else:
        ax3.plot(corr_10.channel, corr_10.counts_corr, marker=".", linestyle="none",
                 label=f"{corr_10.name} (α={corr_10.alpha:.3g})")
        ax3.plot(corr_20.channel, corr_20.counts_corr, marker=".", linestyle="none",
                 label=f"{corr_20.name} (α={corr_20.alpha:.3g})")

    ax3.set_ylabel("counts (bg-subtracted)")
    ax3.set_xlabel("channel")
    ax3.legend()

    fig.suptitle(title)
    plt.tight_layout()
    plt.show()


# -------------------------
# Main (config + run)
# -------------------------

def main() -> None:
    base = Path("../Messdaten/ComptonScattering")

    # --- Configure runs ---
    # Your original script used:
    #   noise: 5 min
    #   signal: 30 min for 10°, 20 min for 20°
    # and then subtracted with factors 6 and 4.
    # Here we compute alpha from times automatically.

    noise_10 = Spectrum(
        name="Noise 10° (5 min)",
        counts=load_tka_counts(base / "Cs137_41B_ring_10_degrees_5min_noise.tka"),
        live_time_s=5 * 60,
    )
    sig_10 = Spectrum(
        name="Ring 10° (30 min)",
        counts=load_tka_counts(base / "Cs137_41B_ring_10_degrees_30min.tka"),
        live_time_s=30 * 60,
    )

    noise_20 = Spectrum(
        name="Noise 20° (5 min)",
        counts=load_tka_counts(base / "Cs137_41B_ring_20_degrees_5min_noise.tka"),
        live_time_s=5 * 60,
    )
    sig_20 = Spectrum(
        name="Ring 20° (20 min)",
        counts=load_tka_counts(base / "Cs137_41B_ring_20_degrees_20min.tka"),
        live_time_s=20 * 60,
    )

    # --- Background subtraction ---
    corr_10 = subtract_background(sig_10, noise_10)  # alpha should be 6
    corr_20 = subtract_background(sig_20, noise_20)  # alpha should be 4

    # --- Plot ---
    plot_comparison(
        noise_10=noise_10,
        noise_20=noise_20,
        sig_10=sig_10,
        sig_20=sig_20,
        corr_10=corr_10,
        corr_20=corr_20,
        show_errorbars_corrected=False,  # set True if you want
    )


if __name__ == "__main__":
    main()