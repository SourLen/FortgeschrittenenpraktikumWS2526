from general_analysis_classes import *
from plotting_functions_ref import *
from gamma_analysis_ref import EnergyCalibration, PeakFit, fit_peak_single, fwhm_energy_from_sigma_channel
from uncertainties import ufloat
import json
#### Load calibration Energy - Channel from energy_calibration.json
with open("./Refactored/energy_calibration.json", "r") as f:
    cal_data = json.load(f)
    energy_cal = EnergyCalibration(
        m=ufloat(cal_data["m"]["value"], cal_data["m"]["uncertainty"]),
        b=ufloat(cal_data["b"]["value"], cal_data["b"]["uncertainty"]),
        cov=np.array(cal_data["cov"])
    )
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

    # Load all spectra into a dictionary
    spectra: dict[str, Spectrum] = {}
    for key, (fname, t) in files.items():
        counts = load_tka_counts(base / fname)
        spectra[key] = Spectrum(name=key, counts=counts, live_time_s=t)

    # Background subtraction for each ring
    ring_corr: dict[str, CorrectedSpectrum] = {}
    for i in [10, 20, 30, 40, 50]:
        signal = spectra[f"ring_{i}"]
        noise = spectra[f"ring_{i}_noise"]
        ring_corr[f"ring_{i}_corr"] = subtract_background(signal, noise)

    ### Plots of all corrected spectra
    fig, ax = plt.subplots(5, 1, figsize=(10, 15))
    for i, angle in enumerate([10, 20, 30, 40, 50]):
        corr_spec = ring_corr[f"ring_{angle}_corr"]
        ax[i].errorbar(corr_spec.channel, corr_spec.counts_corr, yerr=corr_spec.sigma_corr, fmt='o', markersize=2)
        ax[i].set_title(f"Ring at {angle}° (bg-subtracted)")
        ax[i].set_ylabel("Counts")

    plt.xlabel("Channel")
    plt.tight_layout()
    plt.close()
    
    ### Fits on the energy peaks
    ### Peaks: (spectrum, label, (range_lo, range_hi))
    PEAKS = [
            (ring_corr["ring_10_corr"], "Ring 10°", (430, 580)),
            (ring_corr["ring_20_corr"], "Ring 20°", (410, 550)),
            (ring_corr["ring_30_corr"], "Ring 30°", (380, 465)),
            (ring_corr["ring_40_corr"], "Ring 40°", (330, 470)),
            (ring_corr["ring_50_corr"], "Ring 50°", (290, 425)),
        ]
    peakfits: list[PeakFit] = []
    for spec, label, fit_range in PEAKS:
        pf = fit_peak_single(spec, fit_range, label, plot=True)
        peakfits.append(pf)
    ###
    ### Calculate energies from fitted peak positions using the calibration
    angles = []
    energies = []
    energy_errors = []
    for pf in peakfits:
        angle = int(pf.label.split()[1][:-1])  
        energy = energy_cal.energy_of(pf.mu)
        angles.append(angle)
        energies.append(energy.n)
        energy_errors.append(energy.s)
        print(f"Angle: {angle}°, Energy: {energy}")
    
    
    def theoretical_compton_energy(theta_deg: float) -> float:
        """Calculate the theoretical Compton edge energy for a given scattering angle."""
        theta_rad = np.radians(theta_deg)
        E0_keV = 661.8  # Energy of the incident gamma ray in keV (Cs-137)
        m_e_keV = 511.0  # Electron rest mass energy in keV
        return E0_keV / (1 + (E0_keV / m_e_keV) * (1 - np.cos(theta_rad)))
    
    
    
if __name__ == "__main__":
    main()