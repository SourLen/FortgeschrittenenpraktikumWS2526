# gamma_analysis.py
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable, Sequence, Optional
import json

from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties import umath

from activities_ref import SOURCES  # your refactored sources module
from plotting_functions_ref import plot_data_fit_and_pulls
from general_analysis_classes import load_tka_counts, Spectrum, CorrectedSpectrum, subtract_background

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
    ### Fit a single Gaussian + linear background to the specified range of the spectrum.
    spec: CorrectedSpectrum,
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
    s = spec.sigma_corr

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
    spec: CorrectedSpectrum,
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
    cov: np.ndarray 
    def channel_of(self, E_keV: float) -> ufloat:
        return self.m * E_keV + self.b

    def energy_of(self, ch: ufloat) -> ufloat:
        ### Include correlation between m, b
        E = (ch.n - self.b.n) / self.m.n
        # Error propagation with covariance, no correlation b etween ch and m/b since ch is from a different measurement
        dE_dm = -(ch.n - self.b.n) / (self.m.n ** 2)
        dE_db = -1 / self.m.n
        dE_dch = 1 / self.m.n
        var_E = (dE_dm ** 2) * self.cov[0, 0] + (dE_db ** 2) * self.cov[1, 1] + 2 * dE_dm * dE_db * self.cov[0, 1] + (dE_dch ** 2) * ch.s ** 2
        
        return ufloat(E, np.sqrt(var_E))

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
    return EnergyCalibration(m=m, b=b, cov=pcov)

def fwhm_from_sigma_channel(sigma_ch: ufloat) -> ufloat:
    ### FWHM = 2*sqrt(2*ln(2)) * sigma for Gaussian
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
    import matplotlib.pyplot as plt
    base = Path("../Messdaten/Gammaspektroskopie")  # <-- change if needed

    files = {
        "Cs137": ("Cs137_MH851_3min.tka", 180.0),
        "Co60":  ("Co60_LP213_3min.tka", 180.0),
        "Eu152": ("Eu152_MH850_3min.tka", 180.0),
        "Na22":  ("Na22_MH852_3min.tka", 180.0),
        "Noise_5min": ("Noise.tka", 300.0),
        "Noise_D2_5min": ("5min_noise.tka", 300.0),
        "Na22_D2_10min": ("10min_Na22.tka", 600.0),
        "Na22_D2_direct": ("yolo_5min.tka", 300.0),
        "Na22_D2_direct_noise": ("yolo_noise_2min.tka", 120.0),
    }

    spectra: dict[str, Spectrum] = {}
    for key, (fname, t) in files.items():
        counts = load_tka_counts(base / fname)
        spectra[key] = Spectrum(name=key, counts=counts, live_time_s=t)

    # ---- background subtraction ----
    bg = spectra["Noise_5min"]
    cs = subtract_background(spectra["Cs137"], bg)
    co = subtract_background(spectra["Co60"], bg)
    eu = subtract_background(spectra["Eu152"], bg)
    na = subtract_background(spectra["Na22"], bg)

    # D2 special
    bg_d2 = spectra["Noise_D2_5min"]
    na_d2 = subtract_background(spectra["Na22_D2_10min"], bg_d2)
    yolo_na_d2 = spectra["Na22_D2_direct"]
    yolo_bg_d2 = spectra["Na22_D2_direct_noise"]
    yolo_na_d2_corr = subtract_background(yolo_na_d2, bg_d2)
    
    '''
    Additional plot of direct Na22 measurement, showing photopeak but havent yet looked at escape peak- useful for 4.1.5
        
    plt.errorbar(yolo_na_d2_corr.channel, yolo_na_d2_corr.counts_corr, yerr=yolo_na_d2_corr.sigma, fmt=".", label="Na22 D2 direct (BG subtracted)")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.title("Na22 D2 direct measurement (BG subtracted)")
    plt.legend()
    plt.show()
    '''

    # ---- peak definitions ----
    # Each item: (spectrum, label, range_lo, range_hi, literature_energy_keV, I_gamma_fraction, source_key)
    PEAKS = [
        (cs, "Cs-137, 661.7", (450, 550), 661.66, ufloat(0.8500, 0.0020), "Cs-137; MH 851"),
        (co, "Co-60, 1173",   (835, 905), 1173.23, ufloat(0.9985, 0.0003), "Co-60, LP 213"),
        (co, "Co-60, 1332",   (950, 1020),1332.49, ufloat(0.999826, 0.000006), "Co-60, LP 213"),
        (eu, "Eu-152, 121.8", (90, 120),  121.78, ufloat(0.2858, 0.0009), "Eu-152, MH 850"),
        (eu, "Eu-152, 244.7", (170, 215), 244.70, ufloat(0.07580, 0.00030), "Eu-152, MH 850"),
        (eu, "Eu-152, 344.3", (240, 290), 344.28, ufloat(0.265, 0.006), "Eu-152, MH 850"),
        (eu, "Eu-152, 778.9", (545, 620), 778.90, ufloat(0.1294, 0.0015), "Eu-152, MH 850"),
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
    ### Save calibration result to JSON for later use
    cal_dict = {
        "m": {"value": cal.m.n, "uncertainty": cal.m.s},
        "b": {"value": cal.b.n, "uncertainty": cal.b.s},
        "cov": cal.cov.tolist(),
    }
    with open("./Refactored/energy_calibration.json", "w") as f:
        json.dump(cal_dict, f, indent=4)


    # ---- resolution ----
    E_fitted = []
    res_points = []
    # Print resolution points
    for pf, E in zip(peakfits, energies):
        ### R = deltaE/E = 2*sqrt(2*ln(2))* sigma_channel / (mu_channel - b)
        ### Ignoriere Korrelationen zwischen Peakfit und Kalibrierigung
        
        k = 2*np.sqrt(2*np.log(2))
        first_sum = (k**2 * pf.sigma.n**2)/(pf.mu.n - cal.b.n)**4 *(cal.b.s ** 2 + pf.mu.s**2)
        sec_sum = k**2 * pf.sigma.s**2 / (pf.mu.n - cal.b.n)**2
        third_sum = -2*k**2*pf.mu.n*cal.cov[0,1]/(pf.mu.n - cal.b.n)**3
        u_R = np.sqrt(first_sum + sec_sum + third_sum)
        R = k*pf.sigma.n / (pf.mu.n - cal.b.n)
        print(f"{pf.label}: R = {R:.4f} ± {u_R:.4f}")
        res_points.append(ufloat(R, u_R))
        E_fitted.append(cal.energy_of(pf.mu))
        
    def energy_resolution_model(E, a, b):
        return np.sqrt(a**2 + (b**2 / E))
    res_values = np.array([r.n for r in res_points])
    E_fitted_values = np.array([E.n for E in E_fitted])
    popt_res, pcov_res = curve_fit(energy_resolution_model, E_fitted_values, res_values, sigma=np.array([r.s for r in res_points]), absolute_sigma=True)
    plot_data_fit_and_pulls(
        x=E_fitted_values,
        y=res_values,
        yerr=np.array([r.s for r in res_points]),
        model=lambda E, pp: energy_resolution_model(E, pp[0], pp[1]),
        popt=popt_res,
        fit_mask=np.ones_like(E_fitted_values, dtype=bool),
        title="Energy Resolution vs Energy",
        xlabel="Energy (keV)",
        ylabel="Energy Resolution R",
    )
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
if __name__ == "__main__":
    main()