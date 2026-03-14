import numpy as np
import matplotlib.pyplot as plt
from praktikum import analyse
from scipy.optimize import curve_fit

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12
})

# Einlesen der Daten
data_Cs = np.loadtxt("01_T2_Daten/01_Energiespektren/07_Cs137_MH851_5min.tka")
data_Na = np.loadtxt("01_T2_Daten/01_Energiespektren/08_Na22_MH852_5min.tka")
data_Eu = np.loadtxt("01_T2_Daten/01_Energiespektren/06_Eu152_MH850_5min.tka")
data_Co = np.loadtxt("01_T2_Daten/01_Energiespektren/09_Co60_Lp213_5min.tka")
data_NOI = np.loadtxt("01_T2_Daten/01_Energiespektren/01_Rausch_5min.tka")
channels = np.arange(len(data_NOI))
# provisorisch
E_chan = 661.66/330.88

# ---------- Fit-Funktionen ----------
def gauss(x, A, mu, sigma, B):
    return A * np.exp(-(x-mu)**2 / (2*sigma**2)) + B

def fit_peak(x, y, xmin, xmax):
    mask = (x >= xmin) & (x <= xmax)
    x_fit = x[mask]
    y_fit = y[mask]

    # Poisson-Fehler
    y_err = np.sqrt(y_fit)
    y_err[y_err == 0] = 1

    # Startwerte
    A0 = np.max(y_fit) - np.min(y_fit)
    mu0 = x_fit[np.argmax(y_fit)]
    sigma0 = 5
    B0 = np.min(y_fit)

    popt, pcov = curve_fit(
        gauss, x_fit, y_fit,
        p0=[A0, mu0, sigma0, B0],
        sigma=y_err,
        absolute_sigma=True,
        maxfev=20000
    )
    perr = np.sqrt(np.diag(pcov))

    residuals = y_fit - gauss(x_fit, *popt)
    chi2 = np.sum((residuals / y_err)**2)
    dof = len(x_fit) - len(popt)
    chi2_red = chi2 / dof if dof > 0 else np.nan

    return x_fit, y_fit, y_err, popt, perr, residuals, chi2_red

def plot_fit_with_residuals(x_fit, y_fit, y_err, popt, residuals, title):
    fig, (ax1, ax2) = plt.subplots(
        2, 1, sharex=True,
        gridspec_kw={'height_ratios': [3, 1]},
        figsize=(6.5, 6.2)
    )

    ax1.errorbar(x_fit, y_fit, yerr=y_err, fmt='.', label="Daten")
    ax1.plot(x_fit, gauss(x_fit, *popt), label="Gauß-Fit")
    ax1.set_ylabel("Counts")
    ax1.set_title(title)
    ax1.grid(alpha=0.3)
    ax1.legend()

    ax2.axhline(0, lw=1)
    ax2.errorbar(x_fit, residuals, yerr=y_err, fmt='.')
    ax2.set_xlabel("Kanal")
    ax2.set_ylabel("Residuen")
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


# ---------- Deine Daten / Kanäle ----------
# data_Cs, data_Na, data_Eu, data_Co, data_NOI = ... (wie bei dir)
# channels = np.arange(len(data_NOI))

# Optional: Noise abziehen (falls du das für Fits willst)
use_noise_subtracted = False
if use_noise_subtracted:
    spectra = {
        "Cs-137": data_Cs - data_NOI,
        "Na-22": data_Na - data_NOI,
        "Eu-152": data_Eu - data_NOI,
        "Co-60": data_Co - data_NOI,
    }
else:
    spectra = {
        "Cs-137": data_Cs,
        "Na-22": data_Na,
        "Eu-152": data_Eu,
        "Co-60": data_Co,
    }

# ---------- Peak-Jobs: (name, xmin, xmax) pro Quelle ----------
peak_jobs = {
    "Cs-137": [(300, 370)],
    "Co-60":  [(540, 610), (610, 690)],
    "Eu-152": [(60, 80), (160, 190), (370,400),(455,500),(510,570),(660,720)],
    # "Na-22": [(...), (...)]   # falls du später noch Bereiche ergänzt
}
# ---------- Peak-energies: (name, energie) pro Quelle ----------
peak_energies = {
    "Cs-137": [661.66],
    "Co-60":  [1173.2, 1332.5],
    "Eu-152": [121.78, 344.28, 778.90, 964.08, 1085.9, 1112.1]}
    # "Na-22": [(...), (...)]

# ---------- Schleife: fitten + plotten + Ergebnisse sammeln ----------
results = []  # Liste von dicts, easy weiterzuverarbeiten/auszugeben

for source_name, intervals in peak_jobs.items():
    y = spectra[source_name]

    for i, (xmin, xmax) in enumerate(intervals, start=1):
        x_fit, y_fit, y_err, popt, perr, residuals, chi2_red = fit_peak(channels, y, xmin, xmax)

        A, mu, sigma, B = popt
        dA, dmu, dsigma, dB = perr

        title = f"{source_name} – Peak {i} (Fitbereich {xmin}-{xmax})"

        print(f"\n{title}")
        print(f"mu    = {mu:.2f} ± {dmu:.2f}")
        print(f"sigma = {sigma:.2f} ± {dsigma:.2f}")
        print(f"chi2_red = {chi2_red:.2f}")

        plot_fit_with_residuals(x_fit, y_fit, y_err, popt, residuals, title)

        results.append({
            "source": source_name,
            "peak_index": i,
            "xmin": xmin, "xmax": xmax,
            "A": A, "dA": dA,
            "mu": mu, "dmu": dmu,
            "sigma": sigma, "dsigma": dsigma,
            "B": B, "dB": dB,
            "chi2_red": chi2_red,
        })

# Erstes Plotten

# 5 Zeilen, 1 Spalte
fig, axs = plt.subplots(5, 1, figsize=(10,12), sharex=True)

axs[0].step(channels, data_Cs, where="mid")
axs[0].set_title("Cs-137")

axs[1].step(channels, data_Na, where="mid")
axs[1].set_title("Na-22")

axs[2].step(channels, data_Eu, where="mid")
axs[2].set_title("Eu-152")

axs[3].step(channels, data_Co, where="mid")
axs[3].set_title("Co-60")

axs[4].step(channels, data_NOI, where="mid")
axs[4].set_title("Noise")

# Achsenbeschriftung
for ax in axs:
    ax.set_ylabel("Counts")
    ax.grid(alpha=0.3)

axs[-1].set_xlabel("Kanal")

plt.tight_layout()
plt.show()

# Plotten ohne Noise
# 4 Zeilen, 1 Spalte
fig, axs = plt.subplots(4, 1, figsize=(10,12), sharex=True)

axs[0].step(channels, data_Cs - data_NOI, where="mid")
axs[0].set_title("Cs-137")

axs[1].step(channels, data_Na - data_NOI, where="mid")
axs[1].set_title("Na-22")

axs[2].step(channels, data_Eu - data_NOI, where="mid")
axs[2].set_title("Eu-152")

axs[3].step(channels, data_Co - data_NOI, where="mid")
axs[3].set_title("Co-60")


# Achsenbeschriftung
for ax in axs:
    ax.set_ylabel("Counts ohne Noise")
    ax.grid(alpha=0.3)

axs[-1].set_xlabel("Kanal")

plt.tight_layout()
plt.show()

