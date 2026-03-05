from general_compton_scattering import *

def main() -> None:
    base = Path(__file__).parent.parent.parent / "Messdaten" / "ComptonScattering"/ "KonvGeometrie"
    files = {
        "conv_50_al": ("50/20min_AL_50_deg.tka", 1200.0),
        "conv_50_fe": ("50/20min_FE_50_deg.tka", 1200.0),
        "conv_50_noise": ("50/10min_noise_50_deg.tka", 600.0),
        "conv_60_al": ("60/20min_AL_60_deg.tka", 1200.0),
        "conv_60_fe": ("60/20min_FE_60_deg.tka", 1200.0),
        "conv_60_noise": ("60/10min_noise_60_deg.tka", 600.0),
        "conv_80_al": ("80/20min_AL_80_deg.tka", 1200.0),
        "conv_80_fe": ("80/20min_FE_80_deg.tka", 1200.0),
        "conv_80_noise": ("80/10min_noise_80_deg.tka", 600.0),
        "conv_105_al": ("105/20min_AL_105_deg.tka", 1200.0),
        "conv_105_fe": ("105/20min_FE_105_deg.tka", 1200.0),
        "conv_105_noise": ("105/10min_noise_105_deg.tka", 600.0),
        #"conv_135_al": ("135/20min_AL_135_deg.tka", 1200.0),
        #"conv_135_fe": ("135/20min_FE_135_deg.tka", 1200.0),
        #"conv_135_noise": ("135/10min_noise_135_deg.tka", 300.0),
    }
    spectra: dict[str, Spectrum] = {}
    for key, (fname, t) in files.items():
        counts = load_tka_counts(base / fname)
        spectra[key] = Spectrum(name=key, counts=counts, live_time_s=t)
    conv_corr = {}
    angles = [50, 60, 80, 105] #, 135]
    for angle in angles: #, 135]:
        signal_al = spectra[f"conv_{angle}_al"]
        signal_fe = spectra[f"conv_{angle}_fe"]
        noise = spectra[f"conv_{angle}_noise"]
        conv_corr[f"conv_{angle}_al_corr"] = subtract_background(signal_al, noise)
        conv_corr[f"conv_{angle}_fe_corr"] = subtract_background(signal_fe, noise)
    for angle in angles: #, 135]:
        plot_spectrum_and_noise(spectra[f"conv_{angle}_al"], spectra[f"conv_{angle}_noise"])
        plot_spectrum_and_noise(spectra[f"conv_{angle}_fe"], spectra[f"conv_{angle}_noise"])
        
if __name__ == "__main__":
    main()
