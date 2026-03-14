from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat, UFloat, umath

from general_analysis_classes import Spectrum, CorrectedSpectrum, PeakFit, load_tka_counts, subtract_background
from gamma_analysis_ref import EnergyCalibration, fit_peak_single
import activities_ref as act
import absorbtion_ref as ab
import ring_geometry_ref as rg
import convential_geometry_ref as cg
from scipy.constants import c as speed_of_light_m_per_s
import pandas as pd
# =========================================================
# Constants
# =========================================================

E0_CS_KEV = ufloat(661.657, 0.003)
E0_CS_MEV = E0_CS_KEV / 1000.0
I_GAMMA_CS = ufloat(0.8510, 0.0020)
M_E_KEV = ufloat(511.0, 0.00001)

Z_BY_MATERIAL = {
    "Al": 13,
    "Fe": 26,
}


# =========================================================
# Small data containers
# =========================================================

@dataclass(frozen=True)
class Measurement:
    spectrum: CorrectedSpectrum
    label: str
    fit_range: tuple[int, int]
    identifier: str
    angle: UFloat
    material: str
    geometry: str  # "ring" or "conv"


@dataclass(frozen=True)
class GeometryFactors:
    r: UFloat          # scatterer -> detector distance (cm)
    r0: UFloat         # source -> scatterer distance (cm)
    F_D: UFloat        # detector / collimator area defining solid angle (cm^2)
    x_before: UFloat   # air before scatterer (cm)
    x_inside: UFloat   # path length inside scatterer before/after scattering (cm)
    x_after: UFloat    # air after scatterer (cm)


@dataclass(frozen=True)
class AnalysisResult:
    measurement: Measurement
    peakfit: PeakFit
    energy_keV: UFloat
    rate_sinv: UFloat
    efficiency: UFloat
    eta: UFloat
    n_electrons: UFloat
    dsdo_cm2: UFloat


# =========================================================
# Load persisted fits
# =========================================================

def load_energy_calibration(path: str | Path) -> EnergyCalibration:
    with open(path, "r", encoding="utf-8") as f:
        cal_data = json.load(f)
    return EnergyCalibration(
        m=ufloat(cal_data["m"]["value"], cal_data["m"]["uncertainty"]),
        b=ufloat(cal_data["b"]["value"], cal_data["b"]["uncertainty"]),
        cov=np.array(cal_data["cov"]),
    )


def load_efficiency_fit(path: str | Path) -> tuple[np.ndarray, np.ndarray, float | None]:
    with open(path, "r", encoding="utf-8") as f:
        eff_data = json.load(f)
    popt = np.array(eff_data["popt"], dtype=float)
    pcov = np.array(eff_data["pcov"], dtype=float)
    chi2_red = eff_data.get("chi2_red")
    return popt, pcov, chi2_red


# =========================================================
# Helper physics functions
# =========================================================

def theoretical_compton_energy(theta_deg: float) -> float:
    theta_rad = np.radians(theta_deg)
    return E0_CS_KEV.n / (1.0 + (E0_CS_KEV.n / M_E_KEV.n) * (1.0 - np.cos(theta_rad)))


def klein_nishina_mb_per_sr(theta_deg: float, E0_keV: float = E0_CS_KEV.n) -> float:
    """Klein-Nishina differential cross section in milibarn/sr per electron."""
    alpha_fs = 1 / 137.035999084
    lambda_c_cm = 2.42631023867e-10  # electron Compton wavelength in cm
    theta = np.radians(theta_deg)
    a = E0_keV / 511.0
    rho = 1.0 + a * (1.0 - np.cos(theta))
    pref = alpha_fs**2 * lambda_c_cm**2 / (8.0 * np.pi**2)
    return pref * (1.0 / rho**2) * (rho + 1.0 / rho - np.sin(theta)**2)*10**27  # convert cm^2 to mbarn


def efficiency_from_linear_fit(E_keV: UFloat, popt: np.ndarray, pcov: np.ndarray) -> UFloat:
    """epsilon(E) = a E + b with propagated uncertainty incl. fit covariance and E uncertainty."""
    a, b = popt
    eps_nom = a * E_keV.n + b
    var_eps = (
        (E_keV.n**2) * pcov[0, 0]
        + pcov[1, 1]
        + 2.0 * E_keV.n * pcov[0, 1]
        + (a**2) * (E_keV.s**2)
    )
    var_eps = max(var_eps, 0.0)
    return ufloat(eps_nom, np.sqrt(var_eps))


# =========================================================
# Geometry
# =========================================================

def ring_geometry_factors(angle_deg: float, diameter: UFloat = rg.d) -> GeometryFactors:
    idx = int(round(angle_deg / 10.0) - 1)
    if idx < 0 or idx >= len(rg.a):
        raise ValueError(f"Unsupported ring angle {angle_deg}°")

    # Distances from point-like approximation based on measured a,b,d
    r = umath.sqrt(rg.a[idx] ** 2 + (diameter / 2) ** 2)
    r0 = umath.sqrt(rg.b[idx] ** 2 + (diameter / 2) ** 2)

    # Ring geometry path lengths from ring_geometry_ref / absorbtion_ref convention
    x_after = rg.x1[idx]
    x_before = rg.x2[idx]
    x_inside = rg.xring

    F_D = np.pi * (cg.Dz / 2) ** 2
    return GeometryFactors(r=r, r0=r0, F_D=F_D, x_before=x_before, x_inside=x_inside, x_after=x_after)


def conventional_geometry_factors() -> GeometryFactors:
    r = cg.s1 + cg.s2 
    r0 = cg.rT - cg.s0
    x_before = cg.x_luft1
    x_after = cg.x_luft2
    x_inside = cg.x_material
    F_D = np.pi * (cg.Dk / 2) ** 2
    return GeometryFactors(r=r, r0=r0, F_D=F_D, x_before=x_before, x_inside=x_inside, x_after=x_after)


# =========================================================
# Material / absorption / cross section
# =========================================================

def material_dataset(material: str):
    if material == "Al":
        return ab.al_data
    if material == "Fe":
        return ab.fe_data
    raise ValueError(f"Unknown material '{material}'")


def electron_density_per_cm3(material: str) -> float:
    if material == "Al":
        return ab.atom_density_al * Z_BY_MATERIAL[material]
    if material == "Fe":
        return ab.atom_density_fe * Z_BY_MATERIAL[material]
    raise ValueError(f"Unknown material '{material}'")


def scatterer_volume_cm3(material: str, geometry: str, angle_is_10:bool = False) -> UFloat:
    if geometry == "ring" and not angle_is_10:
        # Torus with square cross section: V = 2*pi*R*a^2, with R = d/2 and a = d2
        return 2.0 * np.pi * (rg.d / 2) * rg.d2**2
    elif geometry == "ring" and angle_is_10:
        # Used the smaller ring
        r10_diam = ufloat(13.5, np.sqrt(2 * 0.1**2 / 12))
        return 2.0 * np.pi * (r10_diam/2) * rg.d2**2
    if geometry == "conv":
        # Cylinder: V = pi * (Ds/2)^2 * h
        return np.pi * (cg.Ds / 2) ** 2 * cg.h
    raise ValueError(f"Unknown geometry '{geometry}'")


def electron_count(material: str, geometry: str, angle_is_10: bool = False) -> UFloat:
    return scatterer_volume_cm3(material, geometry, angle_is_10) * electron_density_per_cm3(material)


def absorption_factor(material: str, E_after_keV: UFloat, geo: GeometryFactors) -> UFloat:
    # absorbtion_ref expects energies in MeV
    E_after_MeV = E_after_keV / 1000.0
    return ab.get_absorption(
        E=E_after_MeV,
        E0=E0_CS_MEV,
        x_before=geo.x_before,
        x_after=geo.x_after,
        x_inside=geo.x_inside,
        material=material_dataset(material),
    )


def differential_cross_section(
    rate_sinv: UFloat,
    activity_kBq: UFloat,
    I_gamma: UFloat,
    efficiency: UFloat,
    eta: UFloat,
    n_electrons: UFloat,
    geo: GeometryFactors,
) -> UFloat:
    activity_Bq = activity_kBq * 1000.0
    solid_angle = geo.F_D / geo.r**2
    flux_at_scatterer = activity_Bq * I_gamma / (4.0 * np.pi * geo.r0**2)
    denom = flux_at_scatterer * eta * efficiency * n_electrons * solid_angle
    return rate_sinv / denom


# =========================================================
# Data loading / measurements
# =========================================================

def load_ring_spectra(base: Path) -> dict[str, CorrectedSpectrum]:
    files = {
        "ring_10": ("ring_10_degrees_20min.tka", 1200.0),
        "ring_20": ("ring_20_degrees_20min.tka", 1200.0),
        "ring_30": ("ring_30_degrees_20min.tka", 1200.0),
        "ring_40": ("ring_40_degrees_20min.tka", 1200.0),
        "ring_50": ("ring_50_degrees_20min.tka", 1200.0),
        "ring_10_noise": ("ring_10_degrees_5min_noise.tka", 300.0),
        "ring_20_noise": ("ring_20_degrees_5min_noise.tka", 300.0),
        "ring_30_noise": ("ring_30_degrees_5min_noise.tka", 300.0),
        "ring_40_noise": ("ring_40_degrees_5min_noise.tka", 300.0),
        "ring_50_noise": ("ring_50_degrees_5min_noise.tka", 300.0),
    }

    spectra: dict[str, Spectrum] = {}
    for key, (fname, t) in files.items():
        spectra[key] = Spectrum(name=key, counts=load_tka_counts(base / fname), live_time_s=t)

    corrected: dict[str, CorrectedSpectrum] = {}
    for angle in [10, 20, 30, 40, 50]:
        corrected[f"ring_{angle}_corr"] = subtract_background(spectra[f"ring_{angle}"], spectra[f"ring_{angle}_noise"])
    return corrected


def load_conv_spectra(base: Path) -> dict[str, CorrectedSpectrum]:
    files = {
        "conv_50_al": ("50/20min_Al_50_deg.tka", 1200.0),
        "conv_50_fe": ("50/20min_Fe_50_deg.tka", 1200.0),
        "conv_50_noise": ("50/10min_noise_50_deg.tka", 600.0),
        "conv_60_al": ("60/20min_Al_60_deg.tka", 1200.0),
        "conv_60_fe": ("60/20min_Fe_60_deg.tka", 1200.0),
        "conv_60_noise": ("60/10min_noise_60_deg.tka", 600.0),
        "conv_80_al": ("80/20min_Al_80_deg.tka", 1200.0),
        "conv_80_fe": ("80/20min_Fe_80_deg.tka", 1200.0),
        "conv_80_noise": ("80/10min_noise_80_deg.tka", 600.0),
        "conv_105_al": ("105/20min_Al_105_deg.tka", 1200.0),
        "conv_105_fe": ("105/20min_Fe_105_deg.tka", 1200.0),
        "conv_105_noise": ("105/10min_noise_105_deg.tka", 600.0),
        "conv_135_al": ("135/20min_Al_135_deg.tka", 1200.0),
        "conv_135_fe": ("135/20min_Fe_135_deg.tka", 1200.0),
        "conv_135_noise": ("135/10min_noise_135_deg.tka", 600.0),
    }

    spectra: dict[str, Spectrum] = {}
    for key, (fname, t) in files.items():
        spectra[key] = Spectrum(name=key, counts=load_tka_counts(base / fname), live_time_s=t)

    corrected: dict[str, CorrectedSpectrum] = {}
    for angle in [50, 60, 80, 105, 135]:
        noise = spectra[f"conv_{angle}_noise"]
        corrected[f"conv_{angle}_al_corr"] = subtract_background(spectra[f"conv_{angle}_al"], noise)
        corrected[f"conv_{angle}_fe_corr"] = subtract_background(spectra[f"conv_{angle}_fe"], noise)
    return corrected


def build_measurements(ring_corr: dict[str, CorrectedSpectrum], conv_corr: dict[str, CorrectedSpectrum]) -> list[Measurement]:
    unc_angles_conv = 5 / np.sqrt(12)
    return [
        Measurement(ring_corr["ring_10_corr"], "Ring 10°", (445, 550), "r10", rg.theta[0], "Al", "ring"),
        Measurement(ring_corr["ring_20_corr"], "Ring 20°", (420, 550), "r20", rg.theta[1], "Al", "ring"),
        Measurement(ring_corr["ring_30_corr"], "Ring 30°", (380, 465), "r30", rg.theta[2], "Al", "ring"),
        Measurement(ring_corr["ring_40_corr"], "Ring 40°", (330, 470), "r40", rg.theta[3], "Al", "ring"),
        Measurement(ring_corr["ring_50_corr"], "Ring 50°", (290, 425), "r50", rg.theta[4], "Al", "ring"),
        Measurement(conv_corr["conv_50_al_corr"], "Conv 50° Al", (300, 420), "c50al", ufloat(50, unc_angles_conv), "Al", "conv"),
        Measurement(conv_corr["conv_50_fe_corr"], "Conv 50° Fe", (300, 420), "c50fe", ufloat(50, unc_angles_conv), "Fe", "conv"),
        Measurement(conv_corr["conv_60_al_corr"], "Conv 60° Al", (275, 360), "c60al", ufloat(60, unc_angles_conv), "Al", "conv"),
        Measurement(conv_corr["conv_60_fe_corr"], "Conv 60° Fe", (270, 365), "c60fe", ufloat(60, unc_angles_conv), "Fe", "conv"),
        Measurement(conv_corr["conv_80_al_corr"], "Conv 80° Al", (220, 290), "c80al", ufloat(80, unc_angles_conv), "Al", "conv"),
        Measurement(conv_corr["conv_80_fe_corr"], "Conv 80° Fe", (220, 290), "c80fe", ufloat(80, unc_angles_conv), "Fe", "conv"),
        Measurement(conv_corr["conv_105_al_corr"], "Conv 105° Al", (175, 230), "c105al", ufloat(105, unc_angles_conv), "Al", "conv"),
        Measurement(conv_corr["conv_105_fe_corr"], "Conv 105° Fe", (175, 240), "c105fe", ufloat(105, unc_angles_conv), "Fe", "conv"),
        Measurement(conv_corr["conv_135_al_corr"], "Conv 135° Al", (150, 195), "c135al", ufloat(135, unc_angles_conv), "Al", "conv"),
        Measurement(conv_corr["conv_135_fe_corr"], "Conv 135° Fe", (145, 190), "c135fe", ufloat(135, unc_angles_conv), "Fe", "conv"),
    ]


# =========================================================
# Plotting
# =========================================================

def plot_energy_vs_angle(results: list[AnalysisResult]) -> None:
    ring = [r for r in results if r.measurement.geometry == "ring"]
    conv_al = [r for r in results if r.measurement.geometry == "conv" and r.measurement.material == "Al"]
    conv_fe = [r for r in results if r.measurement.geometry == "conv" and r.measurement.material == "Fe"]

    plt.figure(figsize=(8, 5))
    plt.errorbar([r.measurement.angle.n for r in ring], [r.energy_keV.n for r in ring],
                 xerr=[r.measurement.angle.s for r in ring], yerr=[r.energy_keV.s for r in ring],
                 fmt='o', label='Ring')
    plt.errorbar([r.measurement.angle.n for r in conv_al], [r.energy_keV.n for r in conv_al],
                 xerr=[r.measurement.angle.s for r in conv_al], yerr=[r.energy_keV.s for r in conv_al],
                 fmt='o', label='Konventionell, Al')
    plt.errorbar([r.measurement.angle.n for r in conv_fe], [r.energy_keV.n for r in conv_fe],
                 xerr=[r.measurement.angle.s for r in conv_fe], yerr=[r.energy_keV.s for r in conv_fe],
                 fmt='o', label='Konventionell, Fe')

    a_range = np.linspace(0, 140, 400)
    plt.plot(a_range, [theoretical_compton_energy(a) for a in a_range], color='red', label='Theorie')
    plt.xlabel('Streuwinkel (Grad)')
    plt.ylabel('Energie des gestreuten Photons (keV)')
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_dsdo_vs_angle(results: list[AnalysisResult]) -> None:
    ring = [r for r in results if r.measurement.geometry == "ring"]
    conv_al = [r for r in results if r.measurement.geometry == "conv" and r.measurement.material == "Al"]
    conv_fe = [r for r in results if r.measurement.geometry == "conv" and r.measurement.material == "Fe"]
    ### Convert dsdo from cm^2/sr to mbarn/sr for plotting
    cm2mb = 1e27  # 1 cm^2 = 10^27 mbarn
    plt.figure(figsize=(8, 5))
    plt.errorbar([r.measurement.angle.n for r in ring], [r.dsdo_cm2.n*cm2mb for r in ring],
                 xerr=[r.measurement.angle.s for r in ring], yerr=[r.dsdo_cm2.s*cm2mb for r in ring],
                 fmt='o', label='Ring')
    plt.errorbar([r.measurement.angle.n for r in conv_al], [r.dsdo_cm2.n*cm2mb for r in conv_al],
                 xerr=[r.measurement.angle.s for r in conv_al], yerr=[r.dsdo_cm2.s*cm2mb for r in conv_al],
                 fmt='o', label='Konventionell, Al')
    plt.errorbar([r.measurement.angle.n for r in conv_fe], [r.dsdo_cm2.n*cm2mb for r in conv_fe],
                 xerr=[r.measurement.angle.s for r in conv_fe], yerr=[r.dsdo_cm2.s*cm2mb for r in conv_fe],
                 fmt='o', label='Konventionell, Fe')

    a_range = np.linspace(1, 140, 400)
    plt.plot(a_range, [klein_nishina_mb_per_sr(a) for a in a_range], color='red', label='Klein-Nishina')
    plt.xlabel('Streuwinkel (Grad)')
    plt.ylabel(r'Differenzieller Wirkungsquerschnitt $d\sigma/d\Omega$ (mbarn/sr)')
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_inverse_energy_diff_vs_one_minus_cos(results: list[AnalysisResult]) -> None:
    ring = [r for r in results if r.measurement.geometry == "ring"]
    conv_al = [r for r in results if r.measurement.geometry == "conv" and r.measurement.material == "Al"]
    conv_fe = [r for r in results if r.measurement.geometry == "conv" and r.measurement.material == "Fe"]

    plt.figure(figsize=(8, 5))
    for r in ring:
        x = 1 - np.cos(np.radians(r.measurement.angle.n))
        y = 1 / r.energy_keV - 1 / E0_CS_KEV
        xerr = np.sin(np.radians(r.measurement.angle.n)) * np.radians(r.measurement.angle.s)
        plt.errorbar(x, y.n, xerr=xerr, yerr=y.s, fmt='o', label=f'Ring {r.measurement.angle.n:.0f}°')
    for r in conv_al:
        x = 1 - np.cos(np.radians(r.measurement.angle.n))
        y = 1 / r.energy_keV - 1 / E0_CS_KEV
        xerr = np.sin(np.radians(r.measurement.angle.n)) * np.radians(r.measurement.angle.s)
        plt.errorbar(x, y.n, xerr=xerr, yerr=y.s, fmt='o', label=f'Conv Al {r.measurement.angle.n:.0f}°')
    for r in conv_fe:
        x = 1 - np.cos(np.radians(r.measurement.angle.n))
        y = 1 / r.energy_keV - 1 / E0_CS_KEV
        xerr = np.sin(np.radians(r.measurement.angle.n)) * np.radians(r.measurement.angle.s)
        plt.errorbar(x, y.n, xerr=xerr, yerr=y.s, fmt='o', label=f'Conv Fe {r.measurement.angle.n:.0f}°')

    ### Calculate electronmass through the steepnes of a linear fit onto the data
    ### TODO: properly propagate uncertainties from both x and y into the fit, currently only y is considered
    
    def model(x, m, b):
        return m * x + b
    from scipy.optimize import curve_fit
    x_data = []
    y_data = []
    for r in results:
        x = 1 - np.cos(np.radians(r.measurement.angle.n))
        y = 1 / r.energy_keV - 1 / E0_CS_KEV
        x_data.append(x)
        y_data.append(y.n)
    popt, pcov = curve_fit(model, x_data, y_data)
    m_fit = popt[0]
    m_fit_err = np.sqrt(pcov[0, 0])
    b_fit = popt[1]
    b_fit_err = np.sqrt(pcov[1, 1])
    plt.plot(x_data, model(np.array(x_data), m_fit, b_fit), color='red', label=f'Fit: m={m_fit:.3e}±{m_fit_err:.3e} 1/keV, b={b_fit:.3e}±{b_fit_err:.3e} 1/keV')
    plt.xlabel(r'$1 - \cos(\theta)$')
    plt.ylabel(r'$1/E - 1/E_0$ (1/keV)')
    plt.legend()
    plt.tight_layout()
    plt.show()
    m_e_times_c2_keV = 1 / ufloat(m_fit, m_fit_err)
    print(f"Steigung m = {m_fit:.3e} ± {m_fit_err:.3e} 1/keV bzw Elektronenmasse m_e*c^2 = {m_e_times_c2_keV:.1f} keV")
    
# =========================================================
# Main analysis
# =========================================================
data_to_save = {
    "angle_deg": [],
    "material": [],
    "geometry": [],
    "r_0^2": [],
    "r^2": [],
    "rate_sinv": [],
    "activity_kBq": [],
    "I_gamma": [],
    "efficiency": [],
    "eta": [],
    "n_electrons": [],
    "F_D_cm2": [],
    "dsdo_cm2": [],
}
def main() -> None:
    script_dir = Path(__file__).resolve().parent
    energy_cal = load_energy_calibration(script_dir  / 'energy_calibration.json')
    eff_popt, eff_pcov, chi2_red = load_efficiency_fit(script_dir  / 'efficiency_fit_results.json')

    base_ring = script_dir.parent.parent / 'Messdaten' / 'ComptonScattering' / 'Ringgeometrie'
    base_conv = script_dir.parent.parent / 'Messdaten' / 'ComptonScattering' / 'KonvGeometrie'

    ring_corr = load_ring_spectra(base_ring)
    conv_corr = load_conv_spectra(base_conv)
    measurements = build_measurements(ring_corr, conv_corr)

    activity = act.SOURCES['Cs-137_1'].activity_on()
    conv_geo = conventional_geometry_factors()

    print('Loaded efficiency fit:', eff_popt, 'chi2_red =', chi2_red)
    print('Activity on measurement date:', activity, 'kBq')
    print()

    results: list[AnalysisResult] = []
    ind = 0
    for meas in measurements:
        pf = fit_peak_single(
            meas.spectrum,
            meas.fit_range,
            label=meas.label,
            plot=False,
            conv=(meas.geometry == 'conv'),
        )

        energy_keV = energy_cal.energy_of(pf.mu)
        rate_sinv = pf.area / meas.spectrum.live_time_s
        efficiency = efficiency_from_linear_fit(energy_keV, eff_popt, eff_pcov)

        if meas.geometry == "ring" and meas.identifier == "r10":
            r10_diam = ufloat(13.5, np.sqrt(2 * 0.1**2 / 12))
            geo = ring_geometry_factors(meas.angle.n, diameter = r10_diam)
        elif meas.geometry == "ring":
            geo = ring_geometry_factors(meas.angle.n)
        elif meas.geometry == "conv":
            geo = conv_geo
        else:
            raise ValueError(f"Unknown geometry '{meas.geometry}' for measurement '{meas.identifier}'")
        eta = absorption_factor(meas.material, energy_keV, geo)
        n_e = electron_count(meas.material, meas.geometry, angle_is_10=(meas.geometry == "ring" and meas.identifier == "r10"))
        dsdo = differential_cross_section(rate_sinv, activity, I_GAMMA_CS, efficiency, eta, n_e, geo)
        
        results.append(AnalysisResult(
            measurement=meas,
            peakfit=pf,
            energy_keV=energy_keV,
            rate_sinv=rate_sinv,
            efficiency=efficiency,
            eta=eta,
            n_electrons=n_e,
            dsdo_cm2=dsdo,
        ))

        #print(f'{meas.label:14s}E={energy_keV} keV  activity = {activity}   rate={rate_sinv}  eps={efficiency}  eta={eta}  N_e={n_e}   solid_angle = {geo.F_D / geo.r**2} dsdo={dsdo} cm^2/sr')
        ### Save: angle, material,r_0^2, r^2, rate, activity, I_gamma, eps, eta, N_electron, F_D, dsdo into a csv file
        data_to_save["angle_deg"].append(meas.angle.n)
        data_to_save["material"].append(meas.material)
        data_to_save["geometry"].append(meas.geometry)
        data_to_save["r_0^2"].append(geo.r0.n**2)
        data_to_save["r^2"].append(geo.r.n**2)
        data_to_save["rate_sinv"].append(rate_sinv.n)
        data_to_save["activity_kBq"].append(activity.n)
        data_to_save["I_gamma"].append(I_GAMMA_CS.n)
        data_to_save["efficiency"].append(efficiency.n)
        data_to_save["eta"].append(eta.n)
        data_to_save["n_electrons"].append(n_e.n)
        data_to_save["F_D_cm2"].append(geo.F_D.n)
        data_to_save["dsdo_cm2"].append(dsdo.n)
        
        ind += 1
        print(f"angle: {meas.angle.n}, rate:{rate_sinv}")
    pd.DataFrame(data_to_save).to_csv(script_dir / "compton_scattering_results.csv", sep = ";", decimal = ",", index=False, encoding="utf-8")
    plot_energy_vs_angle(results)
    plot_dsdo_vs_angle(results)
    plot_inverse_energy_diff_vs_one_minus_cos(results)

if __name__ == '__main__':
    main()
