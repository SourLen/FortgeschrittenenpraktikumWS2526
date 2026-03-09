from general_analysis_classes import *
from plotting_functions_ref import *
from gamma_analysis_ref import EnergyCalibration, PeakFit, fit_peak_single, fwhm_energy_from_sigma_channel
from uncertainties import ufloat, UFloat, umath
import json
import activities_ref as act
import absorbtion_ref as ab
import ring_geometry_ref as rg
import convential_geometry_ref as cg

#### Load calibration Energy - Channel from energy_calibration.json
with open("./Refactored/energy_calibration.json", "r") as f:
    cal_data = json.load(f)
    energy_cal = EnergyCalibration(
        m=ufloat(cal_data["m"]["value"], cal_data["m"]["uncertainty"]),
        b=ufloat(cal_data["b"]["value"], cal_data["b"]["uncertainty"]),
        cov=np.array(cal_data["cov"])
    )
    
def comparison_plot(spec1: CorrectedSpectrum, spec2: CorrectedSpectrum, spec3: CorrectedSpectrum, label1: str, label2: str, label3: str) -> None:
    plt.figure(figsize=(10, 5))
    plt.errorbar(spec1.channel, spec1.counts_corr, yerr=spec1.sigma_corr, fmt='o', markersize=2, label=label1)
    plt.errorbar(spec2.channel, spec2.counts_corr, yerr=spec2.sigma_corr, fmt='o', markersize=2, label=label2)
    plt.errorbar(spec3.channel, spec3.counts_corr, yerr=spec3.sigma_corr, fmt='o', markersize=2, label=label3)
    plt.xlabel("Channel")
    plt.ylabel("Counts (bg-subtracted)")
    plt.title(f"Comparison of {label1} and {label2}")
    plt.legend()
    plt.tight_layout()
    #plt.show()
    
def get_r_r0_ring(angle_deg: float) -> tuple[UFloat, UFloat]:
        idx = round(angle_deg)/10 - 1
        r = umath.sqrt(rg.a[int(idx)]**2 + (rg.d/2)**2)
        r0 = umath.sqrt(rg.b[int(idx)]**2 + (rg.d/2)**2)
        return r, r0
def get_eta(material: str, energy_after: UFloat, x_before:UFloat, x_after: UFloat, x_inside: UFloat) -> UFloat:
        if material == "Al":
            mat_data = ab.al_data
        elif material == "Fe":
            mat_data = ab.fe_data
        else:
            raise ValueError("Invalid material")
        return ab.get_absorption(E=energy_after, E0=ufloat(0.661657, 0.000003), x_before=x_before, x_after=x_after, x_inside=x_inside, material=mat_data)
def find_n_e(material: str, is_ring: bool):
        if is_ring:
            volume = 2*np.pi*rg.d2**2*(rg.d / 2)
        else:
            volume = (cg.Ds / 2)**2 * np.pi * cg.h
        if material == "Al":
            N_e = volume * ab.atom_density_al
        elif material == "Fe":
            N_e = volume * ab.atom_density_fe
        else:
            raise ValueError("Invalid material")
        return N_e
     
def diff_cross_section(theta_deg: float, material: str, is_ring: bool, 
                       energy_peak_fitted: UFloat, peak_id: str, kollimDM: UFloat, 
                       detectDM: UFloat, r_conv: UFloat, r0_conv: UFloat, activity: UFloat, I_gamma: UFloat, rates: list[UFloat], eff_points: list[tuple[UFloat, UFloat]],
                       PEAKS: list[tuple[float, float, float, str]]) -> UFloat:
        if is_ring:
            r, r0 = get_r_r0_ring(theta_deg)
            F_D = np.pi*(detectDM/2)**2
            eta = get_eta(material, energy_peak_fitted/1000, x_before=rg.x2[int(theta_deg/10)-1], x_after=rg.x1[int(theta_deg/10)-1], x_inside=rg.xring)
        else:
            r, r0 = r_conv, r0_conv
            F_D = np.pi*(kollimDM/2)**2
            eta = get_eta(material, energy_peak_fitted/1000, x_before=cg.x_luft1, x_after=cg.x_luft2, x_inside=cg.x_material)
        idx = None
        for i in range(len(PEAKS)):
            if PEAKS[i][3] == peak_id:
                idx = i
                break
        if idx is None:
            raise ValueError("Invalid peak_id")
        rate = rates[idx]
        efficiency = eff_points[idx][1]
        print(f"Testeffizienz: {efficiency}")
        N_e = find_n_e(material, is_ring)
        d_sigma = rate / (activity *1e3* I_gamma * (F_D / (4 * np.pi * r**2 * r0**2)) * eta * efficiency * N_e)
        return d_sigma
# -------------------------
# Main (config + run)
# -------------------------

def main() -> None:
    # Use absolute path relative to this script's location
    base = Path(__file__).parent.parent.parent / "Messdaten" / "ComptonScattering"/ "Ringgeometrie"
    # Define files and measurement times for signals and noises
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
    base_conv = Path(__file__).parent.parent.parent / "Messdaten" / "ComptonScattering"/ "KonvGeometrie"
    files_conv = {
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
    
    # Load all spectra into a dictionary
    spectra_ring: dict[str, Spectrum] = {}
    spectr_conv:dict[str, Spectrum] = {}
    for key, (fname, t) in files.items():
        counts = load_tka_counts(base / fname)
        spectra_ring[key] = Spectrum(name=key, counts=counts, live_time_s=t)
    for key, (fname, t) in files_conv.items():
        counts = load_tka_counts(base_conv / fname)
        spectr_conv[key] = Spectrum(name=key, counts=counts, live_time_s=t)
    
    # Background subtraction for each ring and conventional spectrum
    ring_corr: dict[str, CorrectedSpectrum] = {}
    conv_corr: dict[str, CorrectedSpectrum] = {}
    for i in [10, 20, 30, 40, 50]:
        signal = spectra_ring[f"ring_{i}"]
        noise = spectra_ring[f"ring_{i}_noise"]
        ring_corr[f"ring_{i}_corr"] = subtract_background(signal, noise)
    for angle in [50, 60, 80, 105, 135]:
        signal_al = spectr_conv[f"conv_{angle}_al"]
        signal_fe = spectr_conv[f"conv_{angle}_fe"]
        noise = spectr_conv[f"conv_{angle}_noise"]
        conv_corr[f"conv_{angle}_al_corr"] = subtract_background(signal_al, noise)
        conv_corr[f"conv_{angle}_fe_corr"] = subtract_background(signal_fe, noise)
    ### Comparison plot for 50° ring vs conventional geometry
    comparison_plot(ring_corr["ring_50_corr"], conv_corr["conv_50_fe_corr"], conv_corr["conv_50_al_corr"], "Ring 50°", "Conv 50°", "Conv 50° Al")
    
    ### Plot all conventional spectra fe and al each in one subplot, seperated by angles
    fig, ax = plt.subplots(5, 1, figsize=(10, 20), sharex=True)
    angles = [50, 60, 80, 105, 135]
    for i, angle in enumerate(angles):
        ax[i].errorbar(conv_corr[f"conv_{angle}_fe_corr"].channel, conv_corr[f"conv_{angle}_fe_corr"].counts_corr, yerr=conv_corr[f"conv_{angle}_fe_corr"].sigma_corr, fmt='o', markersize=2, label=f"Conv {angle}° Fe")
        ax[i].errorbar(conv_corr[f"conv_{angle}_al_corr"].channel, conv_corr[f"conv_{angle}_al_corr"].counts_corr, yerr=conv_corr[f"conv_{angle}_al_corr"].sigma_corr, fmt='o', markersize=2, label=f"Conv {angle}° Al")
        ax[i].set_ylabel("Counts (bg-subtracted)")
        ax[i].legend()
    ax[-1].set_xlabel("Channel")
    plt.tight_layout()
    #plt.show()
    
    
    
    ### Fits on the energy peaks
    ### Peaks: (spectrum, label, (range_lo, range_hi), identifier)
    PEAKS = [
            (ring_corr["ring_10_corr"], "Ring 10°", (430, 580), "r10"),
            (ring_corr["ring_20_corr"], "Ring 20°", (410, 550), "r20"),
            (ring_corr["ring_30_corr"], "Ring 30°", (380, 465), "r30"),
            (ring_corr["ring_40_corr"], "Ring 40°", (330, 470), "r40"),
            (ring_corr["ring_50_corr"], "Ring 50°", (290, 425), "r50"),
            (conv_corr["conv_50_al_corr"], "Conv 50° Al", (300, 420), "c50al"),
            (conv_corr["conv_50_fe_corr"], "Conv 50° Fe", (300, 420), "c50fe"),
            (conv_corr["conv_60_al_corr"], "Conv 60° Al", (275, 360), "c60al"),
            (conv_corr["conv_60_fe_corr"], "Conv 60° Fe", ( 270, 365), "c60fe"),
            (conv_corr["conv_80_al_corr"], "Conv 80° Al", ( 220, 290), "c80al"),
            (conv_corr["conv_80_fe_corr"], "Conv 80° Fe", ( 225, 290), "c80fe"),
            (conv_corr["conv_105_al_corr"], "Conv 105° Al", (175, 230), "c105al"),
            (conv_corr["conv_105_fe_corr"], "Conv 105° Fe", (175, 240), "c105fe"),
            (conv_corr["conv_135_al_corr"], "Conv 135° Al", (140, 195), "c135al"),
            (conv_corr["conv_135_fe_corr"], "Conv 135° Fe", (145, 190), "c135fe"),
        ]
    peakfits: list[PeakFit] = []
    for spec, label, fit_range, identifier in PEAKS:
        if label.startswith("Conv"):
            pf = fit_peak_single(spec, fit_range, label, plot=False, conv=True)
        else:
            pf = fit_peak_single(spec, fit_range, label, plot=False)
        peakfits.append(pf)
    ###
    ### Calculate energies from fitted peak positions using the calibration
    unc_angles_conv = 5/np.sqrt(12)
    angles = [ufloat(9.99, 0.02), ufloat(20.08, 0.08), ufloat(30.02, 0.2), ufloat(40.01, 0.29), ufloat(49.96, 0.41),
              ufloat(50, unc_angles_conv), ufloat(50, unc_angles_conv),
              ufloat(60, unc_angles_conv), ufloat(60, unc_angles_conv), 
              ufloat(80, unc_angles_conv), ufloat(80, unc_angles_conv),
              ufloat(105, unc_angles_conv), ufloat(105, unc_angles_conv),
              ufloat(135, unc_angles_conv), ufloat(135, unc_angles_conv)]
    energies = []
    energy_errors = []
    for pf in peakfits:
        energy = energy_cal.energy_of(pf.mu)
        energies.append(energy.n)
        energy_errors.append(energy.s)
    
    def theoretical_compton_energy(theta_deg: float) -> float:
        """Calculate the theoretical Compton edge energy for a given scattering angle."""
        theta_rad = np.radians(theta_deg)
        E0_keV = 661.66  # Energy of the incident gamma ray in keV (Cs-137)
        m_e_keV = 511.0  # Electron rest mass energy in keV
        return E0_keV / (1 + (E0_keV / m_e_keV) * (1 - np.cos(theta_rad)))
    plt.errorbar([a.n for a in angles[:5]], energies[:5], xerr = [a.s for a in angles[:5]], yerr=energy_errors[:5], fmt='o', label='Energie (Ring)')
    plt.errorbar([a.n for a in angles[5:][::2]], energies[5:][::2], xerr = [a.s for a in angles[5:][::2]], yerr=energy_errors[5:][::2], fmt='o', label='Energie (Conv, Al)')
    plt.errorbar([a.n for a in angles[6:][::2]], energies[6:][::2], xerr = [a.s for a in angles[6:][::2]], yerr=energy_errors[6:][::2], fmt='o', label='Energie (Conv, Fe)')
    a_range = np.linspace(0, 140, 100)
    plt.plot(a_range, [theoretical_compton_energy(a) for a in a_range], label='Theoretische Vorhersage', color='red')
    plt.xlabel("Streuwinkel (Grad)")
    plt.ylabel("Energie (keV)")
    plt.legend()
    plt.tight_layout()
    #plt.show()
    
    
    ### Rate m per peak (Erst alle Ring, dann Konv abwechselend Al, Fe)
    rates = []
    for (spec, label, fit_range, identifier), pf  in zip(PEAKS, peakfits):
        rates.append(pf.area / spec.live_time_s)
    ### Load efficiency curve from JSON
    with open(".\Refactored\efficiency_fit_results.json", "r") as f:
        eff_data = json.load(f)
    eff_fit_model = lambda E, a, b: a * E + b
    popt = np.array(eff_data["popt"])
    pcov = np.array(eff_data["pcov"])
    ### Calculate efficiency for each peak energy
    eff_points = []
    for energy, energy_error in zip(energies, energy_errors):
        eps_n = eff_fit_model(energy, *popt)
        # Propagate uncertainty from energy to efficiency using the fit parameters
        eps_s = umath.sqrt((popt[0] * energy_error)**2 + (np.sqrt(pcov[0, 0]) * energy)**2 + pcov[1, 1] + 2 * energy * pcov[0, 1])
        eps = ufloat(eps_n, eps_s)
        eff_points.append((energy, eps))
    kollimDM = cg.Dk #cm für konventionelle Geometrie
    detectDM = cg.Dz #cm für Ringgeometrie
    r_conv = cg.s1+cg.s2+cg.rT
    r0_conv = cg.rT-cg.s0
    activity = act.SOURCES["Cs-137_1"].activity_on()
    I_gamma = ufloat(0.8500, 0.0020)
    d_sigmas = []
    for (spec, label, fit_range, identifier), energy_peak_fitted in zip(PEAKS, energies):
        if label.startswith("Ring"):
            is_ring = True
            material = "Al"
        else:
            is_ring = False
            if "Al" in label:
                material = "Al"
            elif "Fe" in label:
                material = "Fe"
            else:
                raise ValueError("Invalid material in label")
        d_sigma = diff_cross_section(theta_deg=angles[PEAKS.index((spec, label, fit_range, identifier))].n, 
                                     material=material, is_ring=is_ring, energy_peak_fitted=ufloat(energy_peak_fitted, 
                                                                                                   energy_errors[PEAKS.index((spec, label, fit_range, identifier))]), 
                                     peak_id=identifier, kollimDM=kollimDM, detectDM=detectDM, r0_conv=r0_conv, r_conv=r_conv,
                                     activity=activity, I_gamma=I_gamma, rates=rates, eff_points=eff_points, PEAKS=PEAKS)
        print(f"{label}: d_sigma = {d_sigma} cm^2")
        d_sigmas.append(d_sigma)
    ### Convert to miliBarn
    d_sigmas = [ds * 1e27 for ds in d_sigmas] # cm^2 to mbarn
    ### Plotting d_sigma vs theta for ring and conv geometries, with error bars
    plt.close()
    plt.figure()
    plt.errorbar([a.n for a in angles[:5]], [ds.n for ds in d_sigmas[:5]], xerr = [a.s for a in angles[:5]], yerr=[ds.s for ds in d_sigmas[:5]], fmt='o', label='Ring')
    plt.errorbar([a.n for a in angles[5:][::2]], [ds.n for ds in d_sigmas[5:][::2]], xerr = [a.s for a in angles[5:][::2]], yerr=[ds.s for ds in d_sigmas[5:][::2]], fmt='o', label='Conv, Al')
    plt.errorbar([a.n for a in angles[6:][::2]], [ds.n for ds in d_sigmas[6:][::2]], xerr = [a.s for a in angles[6:][::2]], yerr=[ds.s for ds in d_sigmas[6:][::2]], fmt='o', label='Conv, Fe')
    plt.xlabel("Streuwinkel (Grad)")
    plt.ylabel("Differentieller Wirkungsquerschnitt (cm^2)")
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    
if __name__ == "__main__":
    main()