# gamma_analysis.py
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence, Optional

from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties import umath

from activities_ref import SOURCES  # your refactored sources module
from plotting_functions_ref import plot_data_fit_and_pulls

# -------------------------
# I/O + preprocessing
# -------------------------

def load_tka(path: Path, skip_header_lines: int = 2) -> np.ndarray:
    """Loads a .tka-like file that contains one count per channel; skips header lines."""
    arr = np.loadtxt(path)
    return np.asarray(arr[skip_header_lines:], dtype=float)

@dataclass(frozen=True)
class Spectrum:
    name: str
    counts: np.ndarray          # raw counts per channel
    live_time_s: float          # measurement time of this spectrum
    channel: np.ndarray         # channel axis (0..N-1)

@dataclass(frozen=True)
class BackgroundSubtracted:
    name: str
    counts_corr: np.ndarray     # N_i - alpha B_i
    sigma: np.ndarray           # sqrt(N_i + alpha^2 B_i)
    alpha: float                # T_s / T_b
    channel: np.ndarray
    live_time_s: float

def subtract_background(
    signal: Spectrum,
    background: Spectrum,
) -> BackgroundSubtracted:
    alpha = signal.live_time_s / background.live_time_s
    counts_corr = signal.counts - alpha * background.counts
    sigma = np.sqrt(np.maximum(signal.counts, 0.0) + (alpha**2) * np.maximum(background.counts, 0.0))
    sigma = np.where(sigma == 0, 1.0, sigma)  # avoid zeros
    return BackgroundSubtracted(
        name=signal.name,
        counts_corr=counts_corr,
        sigma=sigma,
        alpha=alpha,
        channel=signal.channel,
        live_time_s=signal.live_time_s,
    )

# -------------------------
# Peak models + fitting
# -------------------------

def gauss_linbg(x: np.ndarray, A: float, mu: float, sig: float, p0: float, p1: float) -> np.ndarray:
    return A * np.exp(-(x - mu) ** 2 / (2 * sig ** 2)) + (p0 + p1 * x)

def two_gauss_shared_sigma_linbg(
    x: np.ndarray,
    A1: float, mu1: float,
    A2: float, mu2: float,
    sig: float,
    p0: float, p1: float
) -> np.ndarray:
    g1 = A1 * np.exp(-(x - mu1) ** 2 / (2 * sig ** 2))
    g2 = A2 * np.exp(-(x - mu2) ** 2 / (2 * sig ** 2))
    return g1 + g2 + (p0 + p1 * x)

@dataclass(frozen=True)
class PeakFit:
    label: str
    model_name: str
    popt: np.ndarray
    pcov: np.ndarray
    mu: ufloat            # centroid in channel (or effective centroid for multiplet)
    sigma: ufloat         # Gaussian sigma in channel
    area: ufloat          # net Gaussian area in counts (sum for multiplet)
    fit_range: tuple[int, int]

def fit_peak_single(
    spec: BackgroundSubtracted,
    fit_range: tuple[int, int],
    label: str,
    plot: bool = False,
    plot_title: str = "",
    context_pad: int = 40,
    bounds: Optional[tuple[Sequence[float], Sequence[float]]] = None,
) -> PeakFit:
    lo, hi = fit_range
    x = spec.channel
    y = spec.counts_corr
    s = spec.sigma

    fit_mask = (x >= lo) & (x <= hi)

    # initial guesses
    y_win = y[fit_mask]
    x_win = x[fit_mask]
    A0 = float(np.max(y_win))
    mu0 = float(x_win[np.argmax(y_win)])
    sig0 = 10.0
    p00 = float(np.median(y_win))
    p10 = 0.0
    p0 = [A0, mu0, sig0, p00, p10]

    if bounds is None:
        # generic safe bounds
        lower = [0.0, lo, 1.0, -np.inf, -np.inf]
        upper = [np.inf, hi, 80.0, np.inf, np.inf]
        bounds = (lower, upper)

    popt, pcov = curve_fit(
        gauss_linbg,
        x[fit_mask],
        y[fit_mask],
        p0=p0,
        sigma=s[fit_mask],
        absolute_sigma=True,
        maxfev=20000,
        bounds=bounds,
    )

    perr = np.sqrt(np.diag(pcov))
    A, mu, sig, p0bg, p1bg = popt
    uA = ufloat(A, perr[0])
    usig = ufloat(sig, perr[2])

    # Gaussian area in counts/bin convention (bin width = 1 channel)
    area = uA * usig * np.sqrt(2.0 * np.pi)

    fit = PeakFit(
        label=label,
        model_name="gauss+linbg",
        popt=np.asarray(popt),
        pcov=np.asarray(pcov),
        mu=ufloat(mu, perr[1]),
        sigma=usig,
        area=area,
        fit_range=fit_range,
    )

    if plot:
        def model_wrapped(xx: np.ndarray, pp: Sequence[float]) -> np.ndarray:
            return gauss_linbg(xx, *pp)
        vlines = [
            (fit.mu.n, "Peak μ"),
            (fit.mu.n - fit.sigma.n * np.sqrt(2*np.log(2)), "HWHM"),
            (fit.mu.n + fit.sigma.n * np.sqrt(2*np.log(2)), "HWHM"),
        ]
        plot_data_fit_and_pulls(
            x=x,
            y=y,
            yerr=s,
            model=model_wrapped,
            popt=fit.popt,
            fit_mask=fit_mask,
            title=plot_title or f"{spec.name}: {label}",
            xlabel="channel",
            ylabel="counts",
            context_pad=context_pad,
            vlines=vlines,
        )

    return fit

def fit_peak_doublet(
    spec: BackgroundSubtracted,
    fit_range: tuple[int, int],
    label: str,
    mu_guesses: tuple[float, float],
    plot: bool = False,
    plot_title: str = "",
    context_pad: int = 40,
) -> PeakFit:
    """
    Two-gaussian shared sigma + linear bg.
    Returns:
      - area = area1 + area2
      - mu is an intensity/area-weighted centroid in channel (for plotting)
    """
    lo, hi = fit_range
    x = spec.channel
    y = spec.counts_corr
    s = spec.sigma
    fit_mask = (x >= lo) & (x <= hi)

    y_win = y[fit_mask]
    x_win = x[fit_mask]

    A0 = float(np.max(y_win))
    sig0 = 12.0
    p00 = float(np.median(y_win))
    p10 = 0.0

    mu1_0, mu2_0 = mu_guesses
    p0 = [0.6*A0, mu1_0, 0.4*A0, mu2_0, sig0, p00, p10]

    lower = [0.0, lo, 0.0, lo, 1.0, -np.inf, -np.inf]
    upper = [np.inf, hi, np.inf, hi, 100.0, np.inf, np.inf]

    popt, pcov = curve_fit(
        two_gauss_shared_sigma_linbg,
        x[fit_mask],
        y[fit_mask],
        p0=p0,
        sigma=s[fit_mask],
        absolute_sigma=True,
        maxfev=50000,
        bounds=(lower, upper),
    )

    perr = np.sqrt(np.diag(pcov))
    A1, mu1, A2, mu2, sig, p0bg, p1bg = popt

    uA1, uA2 = ufloat(A1, perr[0]), ufloat(A2, perr[2])
    usig = ufloat(sig, perr[4])

    area1 = uA1 * usig * np.sqrt(2*np.pi)
    area2 = uA2 * usig * np.sqrt(2*np.pi)
    area = area1 + area2

    # weighted centroid (by fitted areas)
    mu_eff = (area1 * ufloat(mu1, perr[1]) + area2 * ufloat(mu2, perr[3])) / (area1 + area2)

    fit = PeakFit(
        label=label,
        model_name="2gauss(shared σ)+linbg",
        popt=np.asarray(popt),
        pcov=np.asarray(pcov),
        mu=mu_eff,
        sigma=usig,
        area=area,
        fit_range=fit_range,
    )

    if plot:
        def model_wrapped(xx: np.ndarray, pp: Sequence[float]) -> np.ndarray:
            return two_gauss_shared_sigma_linbg(xx, *pp)

        plot_data_fit_and_pulls(
            x=x,
            y=y,
            yerr=s,
            model=model_wrapped,
            popt=fit.popt,
            fit_mask=fit_mask,
            title=plot_title or f"{spec.name}: {label} (doublet)",
            xlabel="channel",
            ylabel="counts",
            context_pad=context_pad,
            vlines=[(fit.mu.n, "μ_eff")],
        )

    return fit

# -------------------------
# Calibration + resolution
# -------------------------

@dataclass(frozen=True)
class EnergyCalibration:
    # channel = m*E + b
    m: ufloat
    b: ufloat

    def channel_of(self, E_keV: float) -> ufloat:
        return self.m * E_keV + self.b

    def energy_of(self, ch: ufloat) -> ufloat:
        return (ch - self.b) / self.m

    def dEdCh(self) -> ufloat:
        # for linear calibration, derivative is constant
        return 1 / self.m

def fit_linear_calibration(energies_keV: np.ndarray, mus_channel: np.ndarray, mu_unc: np.ndarray) -> EnergyCalibration:
    def lin(E, m, b):
        return m*E + b

    popt, pcov = curve_fit(lin, energies_keV, mus_channel, sigma=mu_unc, absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    m = ufloat(popt[0], perr[0])
    b = ufloat(popt[1], perr[1])
    return EnergyCalibration(m=m, b=b)

def fwhm_from_sigma_channel(sigma_ch: ufloat) -> ufloat:
    return ufloat(2.354820045, 0.0) * sigma_ch

def fwhm_energy_from_sigma_channel(cal: EnergyCalibration, sigma_ch: ufloat) -> ufloat:
    # ΔE = (dE/dch) * Δch = (1/m) * FWHM_ch
    return cal.dEdCh() * fwhm_from_sigma_channel(sigma_ch)

###Plot fit results on each peak

# -------------------------
# Efficiency
# -------------------------

@dataclass(frozen=True)
class CollimatorGeometry:
    radius_cm: ufloat     # hole radius
    d0_cm: ufloat         # distance source -> entrance plane
    L_cm: ufloat          # hole length

    def omega_small_angle(self) -> ufloat:
        # Ω ≈ π a^2 / (d0+L)^2
        return np.pi * self.radius_cm**2 / (self.d0_cm + self.L_cm)**2

def photopeak_efficiency(
    peak_area_counts: ufloat,
    live_time_s: float,
    activity_bq: ufloat,
    I_gamma: ufloat,        # fraction (not percent)
    omega: ufloat,
) -> ufloat:
    rate = peak_area_counts / live_time_s
    return rate * 4*np.pi / (activity_bq * I_gamma * omega)


# -------------------------
# Main analysis
# -------------------------

def main() -> None:
    # ---- configure paths ----
    base = Path("../Messdaten/Gammaspektroskopie")  # <-- change if needed

    files = {
        "Cs137": ("Cs137_MH851_3min.tka", 180.0),
        "Co60":  ("Co60_LP213_3min.tka", 180.0),
        "Eu152": ("Eu152_MH850_3min.tka", 180.0),
        "Na22":  ("Na22_MH852_3min.tka", 180.0),
        "Noise_5min": ("Noise.tka", 300.0),
        "Noise_D2_5min": ("5min_noise.tka", 300.0),
        "Na22_D2_10min": ("10min_Na22.tka", 600.0),
    }

    spectra: dict[str, Spectrum] = {}
    for key, (fname, t) in files.items():
        counts = load_tka(base / fname)
        spectra[key] = Spectrum(name=key, counts=counts, live_time_s=t, channel=np.arange(len(counts)))

    # ---- background subtraction ----
    bg = spectra["Noise_5min"]
    cs = subtract_background(spectra["Cs137"], bg)
    co = subtract_background(spectra["Co60"], bg)
    eu = subtract_background(spectra["Eu152"], bg)
    na = subtract_background(spectra["Na22"], bg)

    # D2 special
    bg_d2 = spectra["Noise_D2_5min"]
    na_d2 = subtract_background(spectra["Na22_D2_10min"], bg_d2)

    # ---- peak definitions ----
    # Each item: (spectrum, label, range_lo, range_hi, literature_energy_keV, I_gamma_fraction, source_key)
    # NOTE: set correct I_gamma fractions based on your tables (percent/100).
    PEAKS = [
        (cs, "Cs-137, 661.7", (450, 550), 661.66, ufloat(0.8500, 0.0020), "Cs-137; MH 851"),
        (co, "Co-60, 1173",   (835, 905), 1173.23, ufloat(0.9985, 0.0003), "Co-60, LP 213"),
        (co, "Co-60, 1332",   (950, 1020),1332.49, ufloat(0.999826, 0.000006), "Co-60, LP 213"),
        (eu, "Eu-152, 121.8", (90, 120),  121.78, ufloat(0.2858, 0.0009), "Eu-152, MH 850"),
        (eu, "Eu-152, 244.7", (170, 215), 244.70, ufloat(0.07580, 0.00030), "Eu-152, MH 850"),
        (eu, "Eu-152, 344.3", (240, 290), 344.28, ufloat(0.265, 0.006), "Eu-152, MH 850"),
        (eu, "Eu-152, 778.9", (545, 620), 778.90, ufloat(0.1294, 0.0015), "Eu-152, MH 850"),
        # The problematic region (often a multiplet). Prefer doublet fit or exclude.
        # (eu, "Eu-152 ~1086/1112", (770, 845), 1085.90, ufloat(0.1021, 0.0004), "Eu-152, MH 850"),
    ]

    # ---- fit peaks ----
    peakfits: list[PeakFit] = []
    energies = []
    mu = []
    mu_unc = []
    sigmas = []

    for spec, label, fr, E, I, source_name in PEAKS:
        pf = fit_peak_single(spec, fr, label=label, plot=False, plot_title = f"{label} keV line")
        peakfits.append(pf)
        energies.append(E)
        mu.append(pf.mu.n)
        mu_unc.append(pf.mu.s)
        sigmas.append(pf.sigma)

    energies = np.asarray(energies, float)
    mu = np.asarray(mu, float)
    mu_unc = np.asarray(mu_unc, float)

    # ---- calibration ----
    cal = fit_linear_calibration(energies, mu, mu_unc)
    print(f"Calibration: channel = ({cal.m})*E + ({cal.b})")

    # ---- resolution ----
    # Print resolution points
    for pf, E in zip(peakfits, energies):
        E_fit = cal.energy_of(pf.mu)
        dE = fwhm_energy_from_sigma_channel(cal, pf.sigma)
        print(f"{pf.label}: E={E_fit} keV, FWHM={dE} keV, res={dE/E_fit}")

    # ---- efficiency geometry ----
    geom = CollimatorGeometry(
        radius_cm=ufloat(1.2375, 0.005/np.sqrt(12)),  # keep your current estimate
        d0_cm=ufloat(12.1, 0.1/np.sqrt(12)),
        L_cm=ufloat(5.0475, 0.01/np.sqrt(12)),
    )
    omega = geom.omega_small_angle()

    import datetime as dt
    measurement_dt = dt.datetime(2026, 3, 3, 14)

    activity_by_name = {s.name: s.activity_on(date=measurement_dt) for s in SOURCES.values()}  # type: ignore

    # ---- efficiencies ----
    eff_points = []
    for (spec, label, fr, E_lit, I, source_name), pf in zip(PEAKS, peakfits):
        A_bq = activity_by_name[source_name] * 1000 
        eps = photopeak_efficiency(
            peak_area_counts=pf.area,
            live_time_s=spec.live_time_s,
            activity_bq=A_bq,
            I_gamma=I,
            omega=omega,
        )
        eff_points.append((E_lit, eps))
        print(f"Efficiency {label}: {eps}")
    # Fit linear model to efficiency over energy
    def lin_eff(E, a, b):
        return a*E + b
    E_fit = np.linspace(min(energies)*0.9, max(energies)*1.1, 100)
    eff_values = np.array([e.n for _, e in eff_points])
    popt, pcov = curve_fit(lin_eff, energies, eff_values)
    chi2 = np.sum((eff_values - lin_eff(energies, *popt))**2/np.array([e.s for _, e in eff_points])**2)
    ndof = len(eff_values) - len(popt)
    chi2_red = chi2 / ndof if ndof > 0 else float("nan")

    # Plot efficiency and linear regression
    import matplotlib.pyplot as plt
    E_plot = [E for E, _ in eff_points]
    eps_n = [e.n for _, e in eff_points]
    eps_s = [e.s for _, e in eff_points]

    plt.errorbar(E_plot, eps_n, yerr=eps_s, fmt="o")
    plt.plot(E_fit, lin_eff(E_fit, *popt), label=r"Linear fit: $\epsilon = aE + b$" + f"\n$\chi^2_{{red}}={chi2_red:.2f}$")
    plt.xlabel("Energy (keV)")
    plt.ylabel(r"Efficiency $\varepsilon$") 
    plt.legend()
    plt.show()
    
    ### Plot all spectra with peak positions and theoretical energies marked, 
    # using the calibration to convert channel to energy on the x-axis.
    ### Make all spectra plots into a single figure with 4 subplots (one for each source), and mark the fitted peak positions and their corresponding energies on each plot.
    # Also exclude exceptional noise in the beginning
    fig, ax = plt.subplots(2, 2, figsize=(15, 10))
    axes = ax.flatten()
    
    for spec in [cs, co, eu, na]:
        axes[0].errorbar(spec.channel, spec.counts_corr, yerr=spec.sigma, fmt=".", label=spec.name)
        # Mark fitted peaks
        for pf in peakfits:
            if pf.label.startswith(spec.name.split(",")[0][:2]):
                
                E_fit = cal.energy_of(pf.mu)
                plt.axvline(pf.mu.n, color="red", linestyle="--", label=f"{pf.label} keV expected")
                plt.text(pf.mu.n+10, max(spec.counts_corr)*1.2, f"{E_fit:.1f} keV", rotation=0, verticalalignment='center', color="red")
        plt.xlabel("Channel")
        plt.ylabel("Noise-subtracted counts")
        plt.title(f"Spectrum: {spec.name}")
        plt.legend()
        plt.show()

if __name__ == "__main__":
    main()