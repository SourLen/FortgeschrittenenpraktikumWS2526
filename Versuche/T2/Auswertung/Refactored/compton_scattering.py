from general_compton_scattering import *


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
    ring_corr = {}
    for i in [10, 20, 30, 40, 50]:
        signal = spectra[f"ring_{i}"]
        noise = spectra[f"ring_{i}_noise"]
        ring_corr[f"ring_{i}_corr"] = subtract_background(signal, noise)

    ### Plot all spectra and their corresponding noises in a single figure
    for i in [10, 20, 30, 40, 50]:
        plot_spectrum_and_noise(spectra[f"ring_{i}"], spectra[f"ring_{i}_noise"])
            

if __name__ == "__main__":
    main()