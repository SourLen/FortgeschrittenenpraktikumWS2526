import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from Auswertung_T2_functions import (
    gauss,
    lin,
    res_model,
    fit_peak_with_noise,
    plot_fit_with_residuals_and_lines,
)

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12
})

# =========================================================
# Daten laden
# =========================================================
data_Cs  = np.loadtxt("01_T2_Daten/01_Energiespektren/07_Cs137_MH851_5min.tka")
data_Na  = np.loadtxt("01_T2_Daten/01_Energiespektren/08_Na22_MH852_5min.tka")
data_Eu  = np.loadtxt("01_T2_Daten/01_Energiespektren/06_Eu152_MH850_5min.tka")
data_Co  = np.loadtxt("01_T2_Daten/01_Energiespektren/09_Co60_Lp213_5min.tka")
data_NOI = np.loadtxt("01_T2_Daten/01_Energiespektren/01_Rausch_5min.tka")

channels = np.arange(len(data_NOI))

# =========================================================
# Rohdaten plotten
# =========================================================
fig, axs = plt.subplots(5, 1, figsize=(10, 12), sharex=True)

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

for ax in axs:
    ax.set_ylabel("Counts")
    ax.grid(alpha=0.3)

axs[-1].set_xlabel("Kanal")
plt.tight_layout()
plt.show()

# =========================================================
# Spektren nach Noise-Abzug
# =========================================================
fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

axs[0].step(channels, data_Cs - data_NOI, where="mid")
axs[0].set_title("Cs-137")

axs[1].step(channels, data_Na - data_NOI, where="mid")
axs[1].set_title("Na-22")

axs[2].step(channels, data_Eu - data_NOI, where="mid")
axs[2].set_title("Eu-152")

axs[3].step(channels, data_Co - data_NOI, where="mid")
axs[3].set_title("Co-60")

for ax in axs:
    ax.set_ylabel("Counts ohne Noise")
    ax.grid(alpha=0.3)

axs[-1].set_xlabel("Kanal")
plt.tight_layout()
plt.show()

# =========================================================
# Datenstrukturen
# =========================================================
spectrum = {
    "Cs-137": data_Cs,
    "Na-22":  data_Na,
    "Eu-152": data_Eu,
    "Co-60":  data_Co,
}

# Fitintervalle für die einzelnen Peaks
peak_jobs = {
    "Cs-137": [(310, 350)],
    "Co-60":  [(550, 600), (620, 670)],
    "Eu-152": [(65, 85), (166, 196), (360, 410), (445, 495)],
    # "Na-22": [(250, 270)]
}

# Zugehörige Literaturenergien in derselben Reihenfolge wie peak_jobs
peak_energies = {
    "Cs-137": [661.66],
    "Co-60":  [1173.2, 1332.5],
    "Eu-152": [121.78, 344.28, 778.90, 964.08],
    # "Na-22": [511.0]
}

activities = {
    "Cs-137": (22310, 90),
    "Na-22":  (10.601, 0.26),
    "Eu-152": (12001, 1.4),
    "Co-60":  (1824.9, 1.1),
}

gamma_rates = {
    "Cs-137": [(0.85, 0.002)],
    "Na-22":  [1.798],  # aktuell im weiteren Code nicht benutzt
    "Eu-152": [(0.2858, 0.0009), (0.265, 0.006), (0.1294, 0.0015), (0.1460, 0.0004)],
    "Co-60":  [(0.9985, 0.0003), (0.999826, 0.000006)],
}

# =========================================================
# Peak-Fits
# =========================================================
results = []

for source_name, intervals in peak_jobs.items():
    y = spectrum[source_name]

    for i, (xmin, xmax) in enumerate(intervals, start=1):
        x_fit, y_fit, y_err, popt, perr, residuals, chi2_red = fit_peak_with_noise(
            channels, y, data_NOI, xmin, xmax
        )

        A, mu, sigma, B = popt
        dA, dmu, dsigma, dB = perr

        # Literaturwert des Peaks passend zum Intervall
        E_theo = peak_energies[source_name][i - 1]

        title = f"{source_name} – Peak {i} (Fitbereich {xmin}-{xmax})"

        print(f"\n{title}")
        print(f"mu    = {mu:.2f} ± {dmu:.2f} (Kanal)")
        print(f"sigma = {sigma:.2f} ± {dsigma:.2f} (Kanal)")
        print(f"chi2_red = {chi2_red:.2f}")

        # Schwarze Linie hier in der Peak-Voruntersuchung:
        # Mittelpunkt des gewählten Fitbereichs
        mu_ref = 0.5 * (xmin + xmax)

        plot_fit_with_residuals_and_lines(
            x_fit=x_fit,
            y_fit=y_fit,
            y_err=y_err,
            popt=popt,
            residuals=residuals,
            title=title,
            x_theo=mu_ref,
            theo_label="Fitfenster-Mitte"
        )

        results.append({
            "source": source_name,
            "peak_index": i,
            "xmin": xmin,
            "xmax": xmax,
            "A": A,
            "dA": dA,
            "mu": mu,
            "dmu": dmu,
            "sigma": sigma,
            "dsigma": dsigma,
            "B": B,
            "dB": dB,
            "chi2_red": chi2_red,
            "E_lit": E_theo,
        })

# =========================================================
# Kalibration: Kanal(mu) als Funktion der Energie
# =========================================================
E_all, mu_all, dmu_all, sig_all, dsig_all = [], [], [], [], []

for source_name, E_list in peak_energies.items():
    peaks = [r for r in results if r["source"] == source_name]

    for k in range(len(peaks)):
        E_all.append(E_list[k])
        mu_all.append(peaks[k]["mu"])
        dmu_all.append(peaks[k]["dmu"])
        sig_all.append(peaks[k]["sigma"])
        dsig_all.append(peaks[k]["dsigma"])

E_all = np.array(E_all, dtype=float)
mu_all = np.array(mu_all, dtype=float)
dmu_all = np.array(dmu_all, dtype=float)
sig_all = np.array(sig_all, dtype=float)
dsig_all = np.array(dsig_all, dtype=float)

dmu_all[dmu_all == 0] = 1.0
dsig_all[dsig_all == 0] = 1.0

# gewichteter linearer Fit: mu = m*E + q
p0 = [
    (mu_all.max() - mu_all.min()) / (E_all.max() - E_all.min()),
    mu_all.mean() - ((mu_all.max() - mu_all.min()) / (E_all.max() - E_all.min())) * E_all.mean()
]

popt, pcov = curve_fit(
    lin,
    E_all,
    mu_all,
    p0=p0,
    sigma=dmu_all,
    absolute_sigma=True,
    maxfev=20000
)

perr = np.sqrt(np.diag(pcov))
m, q = popt
dm, dq = perr

mu_fit = lin(E_all, m, q)
residuals = mu_all - mu_fit

chi2 = np.sum((residuals / dmu_all) ** 2)
dof = len(E_all) - len(popt)
chi2_red = chi2 / dof if dof > 0 else np.nan

print("\n====================")
print("Energie-Kalibration (linear): mu = m*E + q")
print(f"m = {m:.6f} ± {dm:.6f} Kanal/keV")
print(f"q = {q:.3f} ± {dq:.3f} Kanal")
print(f"chi2_red = {chi2_red:.2f}")
print("====================\n")

# =========================================================
# Plot Kalibration + Residuen
# =========================================================
order = np.argsort(E_all)
E_s = E_all[order]
mu_s = mu_all[order]
dmu_s = dmu_all[order]

E_grid = np.linspace(E_s.min() * 0.98, E_s.max() * 1.02, 400)

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True,
    gridspec_kw={"height_ratios": [3, 1]},
    figsize=(7.2, 6.4)
)

ax1.errorbar(E_s, mu_s, yerr=dmu_s, fmt='.', label="Peaks (μ ± dμ)")
ax1.plot(E_grid, lin(E_grid, m, q), label=f"Fit: μ = {m:.6f}·E + {q:.2f}")
ax1.set_ylabel("Kanal μ")
ax1.set_title("Kalibration: Kanal als Funktion der Energie")
ax1.grid(alpha=0.3)
ax1.legend(loc="upper right")

ax2.axhline(0, lw=1)
ax2.errorbar(E_all, residuals, yerr=dmu_all, fmt='.')
ax2.set_xlabel("Energie E (keV)")
ax2.set_ylabel("Residuen (Kanal)")
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# =========================================================
# Halbwertsbreiten / Energieauflösung
# =========================================================
k = 2 * np.sqrt(np.log(4))

# FWHM in Channels
FWHM = k * sig_all
dFWHM = k * dsig_all

# FWHM in keV
FWHM_E = FWHM / m
dFWHM_E = np.sqrt((dFWHM / m) ** 2 + (FWHM * dm / m**2) ** 2)

FWHM_E_s = FWHM_E[order]
dFWHM_E_s = dFWHM_E[order]

plt.figure()
plt.errorbar(E_s, FWHM_E_s, yerr=dFWHM_E_s, fmt='.')
plt.xlabel("Energie (keV)")
plt.ylabel("FWHM (keV)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

R = FWHM_E_s / E_s
dR = dFWHM_E_s / E_s

p0_nl = [0.01, 1.0]

popt_nl, pcov_nl = curve_fit(
    res_model,
    E_s,
    R,
    p0=p0_nl,
    sigma=dR,
    absolute_sigma=True,
    maxfev=20000,
    bounds=(0, np.inf)
)

perr_nl = np.sqrt(np.diag(pcov_nl))
a_nl, b_nl = popt_nl
da_nl, db_nl = perr_nl

print("\n====================")
print("Nichtlinearer Fit: dE/E = sqrt(a^2 + b^2/E)")
print(f"a = {a_nl:.5f} ± {da_nl:.5f}")
print(f"b = {b_nl:.5f} ± {db_nl:.5f}  [sqrt(keV)]")
print("====================\n")

E_grid = np.linspace(E_s.min() * 0.95, E_s.max() * 1.05, 400)
R_fit = res_model(E_s, a_nl, b_nl)
resid = R - R_fit

chi2 = np.sum((resid / dR) ** 2)
dof = len(E_s) - 2
chi2_red = chi2 / dof if dof > 0 else np.nan

print(f"chi2_red (nichtlinear) = {chi2_red:.2f}")

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, figsize=(7.2, 6.4),
    gridspec_kw={"height_ratios": [3, 1]}
)

ax1.errorbar(E_s, R, yerr=dR, fmt='.', linestyle='none', label="Daten: FWHM/E")
ax1.plot(E_grid, res_model(E_grid, a_nl, b_nl), label=f"Fit: a={a_nl:.4f}, b={b_nl:.4f}")
ax1.set_ylabel("dE/E")
ax1.set_title("Energieauflösung: nichtlinearer Fit")
ax1.grid(alpha=0.3)
ax1.legend(loc="upper right")

ax2.axhline(0, lw=1)
ax2.errorbar(E_s, resid, yerr=dR, fmt='.', linestyle='none')
ax2.set_xlabel("Energie E (keV)")
ax2.set_ylabel("Residuen")
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# =========================================================
# Detektoreffizienz
# eps = M / (A * I * Omega_D)
# =========================================================
epsilons = []

r0 = 0.2575
dr0 = 0.001 / np.sqrt(12)

Omega_D = 2.479 * 0.001
dOmega_D = 0.00001

t = 300  # 5 min

for source_name, I in gamma_rates.items():
    peaks = [r for r in results if r["source"] == source_name]

    for k in range(len(peaks)):
        M = np.sqrt(2 * np.pi) * peaks[k]["sigma"] * peaks[k]["A"] / t

        dM = M * np.sqrt(
            (peaks[k]["dA"] / peaks[k]["A"]) ** 2 +
            (peaks[k]["dsigma"] / peaks[k]["sigma"]) ** 2
        )

        Eps = M / (activities[source_name][0] * I[k][0] * Omega_D)

        rel2 = (
            (dM / M) ** 2
            + (2 * dr0 / r0) ** 2
            + (activities[source_name][1] / activities[source_name][0]) ** 2
            + (I[k][1] / I[k][0]) ** 2
            + (dOmega_D / Omega_D) ** 2
        )

        dEps = Eps * np.sqrt(rel2)

        print(source_name)
        print(activities[source_name][0])
        print(I[k][0])

        epsilons.append((source_name, k + 1, Eps, dEps))

E_eps, eps, deps = [], [], []

for source_name, peak_idx, Eps, dEps in epsilons:
    E0 = peak_energies[source_name][peak_idx - 1]
    E_eps.append(E0)
    eps.append(Eps)
    deps.append(dEps)

E_eps = np.array(E_eps, float)
eps = np.array(eps, float)
deps = np.array(deps, float)

ordr = np.argsort(E_eps)

plt.figure()
plt.errorbar(E_eps[ordr], eps[ordr], yerr=deps[ordr], fmt='.', linestyle='none')
plt.xlabel(r"E$_\gamma$ (keV)")
plt.ylabel(r"$\varepsilon$")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()