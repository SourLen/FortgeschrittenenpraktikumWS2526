import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from Auswertung_T2_functions import *
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12
})

# -------------------- Daten laden --------------------

data_Cs = np.loadtxt("01_T2_Daten/01_Energiespektren/07_Cs137_MH851_5min.tka")
data_Na = np.loadtxt("01_T2_Daten/01_Energiespektren/08_Na22_MH852_5min.tka")
data_Eu = np.loadtxt("01_T2_Daten/01_Energiespektren/06_Eu152_MH850_5min.tka")
data_Co = np.loadtxt("01_T2_Daten/01_Energiespektren/09_Co60_Lp213_5min.tka")
data_NOI = np.loadtxt("01_T2_Daten/01_Energiespektren/01_Rausch_5min.tka")


channels = np.arange(len(data_NOI))

# -------------------- Dictionaries --------------------
spectrum = {
    "Cs-137": data_Cs,
    "Na-22": data_Na,
    "Eu-152": data_Eu,
    "Co-60": data_Co,
    }

peak_jobs = {
    "Cs-137": [(310, 350)],
    "Co-60":  [(550, 600), (620, 670)],
    "Eu-152": [(65, 85), (166, 196), (360, 410), (445, 495),],# (510, 570)],
    #"Na-22": [(250,270)]
}

peak_energies = {
    "Cs-137": [661.66],
    "Co-60":  [1173.2, 1332.5],
    "Eu-152": [121.78, 344.28, 778.90, 964.08,]# 1085.9, 1112.1],
    #"Na-22": [511]
}
 
activities = {
    "Cs-137": (22.30e3, 0.09e3),
    "Na-22": (0.10570e3, 0.00026e3),
    "Eu-152": (11.995e3, 0.014e3),
    "Co-60": (1.8222e3, 0.0011e3),
}

gamma_rates = {
    "Cs-137": [(0.85,0.002)] ,
    "Na-22": [1.798],  
    "Eu-152": [(0.2858, 0.0009), (0.265, 0.006), (0.1294, 0.0015), (0.1460,0.0004)],#, (0.1021, 0.0004), (0.1364,0.0004)],
    "Co-60": [(0.9985, 0.0003), (0.999826, 0.000006)],
    }

# -------------------- Erstes Plotten --------------------
# Rohdaten
plot_spectra(
    channels,
    spectrum,
    noise=data_NOI,
    subtract_noise=False,
    ylabel="Counts"
)

# Ohne Noise
plot_spectra(
    channels,
    spectrum,
    noise=data_NOI,
    subtract_noise=True,
    ylabel="Counts ohne Noise"
)


# -------------------- Loop: Peaks fitten, Ergebnisse sammeln --------------------
results = []

for source_name, intervals in peak_jobs.items():
    y = spectrum[source_name]
     
    
    for i, (xmin, xmax) in enumerate(intervals, start=1):
        x_fit, y_fit, y_err, popt, perr, residuals, chi2_red = fit_peak_with_noise(channels, y, data_NOI, xmin, xmax)
        
        A, mu, sigma, B = popt
        dA, dmu, dsigma, dB = perr
        
        title = f"{source_name} – Peak {i} (Fitbereich {xmin}-{xmax}) , $\\mu = {mu:.2f} \\pm {dmu:.2f} , \\chi^2_\\mathrm{{red}} = {chi2_red:.2f}$"

        print(f"\n{title}")
        print(f"mu    = {mu:.2f} ± {dmu:.2f} (Kanal)")
        print(f"sigma = {sigma:.2f} ± {dsigma:.2f} (Kanal)")
        print(f"chi2_red = {chi2_red:.2f}")
        
        plot_fit_with_residuals_and_lines(x_fit, y_fit, y_err, popt, residuals, title)
        
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


# -------------------- Kalibration: Kanal(mu) als Funktion von Energie --------------------

# Arrays (E_all, mu_all, dmu_all) in korrekter Reihenfolge bauen
E_all, mu_all, dmu_all,sig_all, dsig_all = [], [], [], [], []

for source_name, E_list in peak_energies.items():
    peaks = [r for r in results if r["source"] == source_name]

    n = len(peaks)

    for k in range(n):
        E_all.append(E_list[k])
        mu_all.append(peaks[k]["mu"])
        dmu_all.append(peaks[k]["dmu"])
        sig_all.append(peaks[k]["sigma"])
        dsig_all.append(peaks[k]["dsigma"])
        
# konvertierung zu np.arrays
E_all = np.array(E_all, dtype=float)
mu_all = np.array(mu_all, dtype=float)
dmu_all = np.array(dmu_all, dtype=float)
dmu_all[dmu_all == 0] = 1.0  # safety
sig_all = np.array(sig_all, dtype=float)
dsig_all = np.array(dsig_all, dtype=float)
dsig_all[dsig_all == 0] = 1.0  # safety

# Linearer, gewichteter Fit
p0 = [(mu_all.max() - mu_all.min()) / (E_all.max() - E_all.min()),
      mu_all.mean() - ((mu_all.max() - mu_all.min()) / (E_all.max() - E_all.min()))*E_all.mean()]

popt, pcov = curve_fit(
    lin, E_all, mu_all,
    p0=p0,
    sigma=dmu_all,
    absolute_sigma=True,
    maxfev=20000
)

perr = np.sqrt(np.diag(pcov))
m, q = popt
dm, dq = perr

# Residuen und Chi2
mu_fit = lin(E_all, m, q)
residuals = mu_all - mu_fit
chi2 = np.sum((residuals / dmu_all)**2)
dof = len(E_all) - len(popt)
chi2_red = chi2 / dof if dof > 0 else np.nan

print("\n====================")
print("Energie-Kalibration (linear): mu = m*E + q")
print(f"m = {m:.6f} ± {dm:.6f} Kanal/keV")
print(f"q = {q:.3f} ± {dq:.3f} Kanal")
print(f"chi2_red = {chi2_red:.2f}")
print("====================\n")

# -------------------- Plot Kalibration + Residuen --------------------
order = np.argsort(E_all)
E_s = E_all[order]
mu_s = mu_all[order]
dmu_s = dmu_all[order]

E_grid = np.linspace(E_s.min()*0.98, E_s.max()*1.02, 400)

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True,
    gridspec_kw={'height_ratios': [3, 1]},
    figsize=(7.2, 6.4)
)

# Kalibrationsplot
ax1.errorbar(E_s, mu_s, yerr=dmu_s, fmt='.', label="Peaks (μ ± dμ)")
ax1.plot(E_grid, lin(E_grid, m, q), label=f"Fit: μ = {m:.6f}·E + {q:.2f}")
ax1.set_ylabel("Kanal μ")
ax1.set_title("Kalibration: Kanal als Funktion der Energie")
ax1.grid(alpha=0.3)
ax1.legend()

# Residuenplot
ax2.axhline(0, lw=1)
ax2.errorbar(E_all, residuals, yerr=dmu_all, fmt='.')
ax2.set_xlabel("Energie E (keV)")
ax2.set_ylabel("Residuen (Kanal)")
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# -------------------- Zweites Plotten --------------------

# Mit Noise
plot_spectra_with_peaks(
    channels,
    spectrum,
    results,
    peak_energies=peak_energies,
    m=m,
    q=q,
    noise=data_NOI,
    subtract_noise=False,
    ylabel="Counts"
)

# Ohne Noise
plot_spectra_with_peaks(
    channels,
    spectrum,
    results,
    peak_energies=peak_energies,
    m=m,
    q=q,
    noise=data_NOI,
    subtract_noise=True,
    ylabel="Counts ohne Noise"
)

# -------------------- Bestimmen der Halbwertsbreiten --------------------

k = 2 * np.sqrt(np.log(4))
# Als Channel:
FWHM = k * sig_all
dFWHM = k * dsig_all
# Als Energien
FWHM_E = FWHM / m
dFWHM_E = np.sqrt((dFWHM / m)**2 + (FWHM * dm / m**2)**2)
#sortieren fürs plotten
FWHM_E_s = FWHM_E[order]
dFWHM_E_s = dFWHM_E[order]

plt.figure()
plt.errorbar(E_s, FWHM_E_s, yerr=dFWHM_E_s,fmt='.')
plt.xlabel("Energie (keV)")
plt.ylabel("FWHM (keV)")
plt.grid(alpha=0.3)
plt.show()


R = FWHM_E_s / E_s
dR = dFWHM_E_s / E_s

p0_nl = [0.01, 1.0] 

popt_nl, pcov_nl = curve_fit(
    res_model, E_s, R,
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

# Fitkurve + Residuen
E_grid = np.linspace(E_s.min()*0.95, E_s.max()*1.05, 400)
R_fit  = res_model(E_s, a_nl, b_nl)
resid  = R - R_fit

chi2 = np.sum((resid / dR)**2)
dof  = len(E_s) - 2
chi2_red = chi2 / dof if dof > 0 else np.nan
print(f"chi2_red (nichtlinear) = {chi2_red:.2f}")

fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, figsize=(7.2, 6.4),
    gridspec_kw={'height_ratios': [3, 1]}
)

# oben: Daten + Fit
ax1.errorbar(E_s, R, yerr=dR, fmt='.', linestyle='none', label="Daten: FWHM/E")
ax1.plot(E_grid, res_model(E_grid, a_nl, b_nl), label=f"Fit: a={a_nl:.4f}, b={b_nl:.4f}")
ax1.set_ylabel("dE/E")
ax1.set_title("Energieauflösung: nichtlinearer Fit")
ax1.grid(alpha=0.3)
ax1.legend()

# unten: Residuen
ax2.axhline(0, lw=1)
ax2.errorbar(E_s, resid, yerr=dR, fmt='.', linestyle='none')
ax2.set_xlabel("Energie E (keV)")
ax2.set_ylabel("Residuen")
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.show()


# -------------------- Bestimmen der Detektoreffizienz --------------------
# eps = 4 pi * r^2 * F_D / (A*I)
epsilons = []
r0 = 0.2575       # m
dr0 = 0.0042      # m
Fd = 2.03e-3      # m^2
dFd = 0.07e-3     # m^2

Omega_D = Fd / (4*np.pi * r0**2)
dOmega_D = Omega_D * np.sqrt((dFd/Fd)**2 + (2*dr0/r0)**2)

T = {
    "Cs-137": 180,
    "Na-22": 300,
    "Eu-152": 300,
    "Co-60": 300,
}

for source_name, I in gamma_rates.items():
    peaks = [r for r in results if r["source"] == source_name]
    t = T[source_name]

    for k in range(len(peaks)):
        M = np.sqrt(2*np.pi) * peaks[k]["sigma"] * peaks[k]["A"] / t
        dM = M * np.sqrt(
            (peaks[k]["dA"]/peaks[k]["A"])**2 +
            (peaks[k]["dsigma"]/peaks[k]["sigma"])**2
        )

        Eps = M / (activities[source_name][0] * I[k][0] * Omega_D)
        rel2 = (
            (dM/M)**2 +
            (activities[source_name][1]/activities[source_name][0])**2 +
            (I[k][1]/I[k][0])**2 +
            (dOmega_D/Omega_D)**2
        )
        dEps = Eps * np.sqrt(rel2)

        epsilons.append((source_name, k+1, Eps, dEps))

E_eps, eps, deps, labels = [], [], [], []

for (source_name, peak_idx, Eps, dEps) in epsilons:
    E0 = peak_energies[source_name][peak_idx - 1]
    E_eps.append(E0)
    eps.append(Eps)
    deps.append(dEps)
    labels.append(f"{source_name} Peak {peak_idx}\n(E={E0:.0f} keV)")

E_eps = np.array(E_eps, float)
eps = np.array(eps, float)
deps = np.array(deps, float)

# ----- Fit an Effizienz/Energie -----

popt, pcov = curve_fit(
    lin,
    E_eps,
    eps,
    sigma=deps,
    absolute_sigma=True
)

a, b = popt
da, db = np.sqrt(np.diag(pcov))

print("\nLinearer Effizienz-Fit")
print(f"a = {a:.3e} ± {da:.3e} 1/keV")
print(f"b = {b:.3e} ± {db:.3e}")

eps_fit = lin(E_eps, a, b)
residuals = eps - eps_fit

chi2 = np.sum((residuals / deps)**2)
dof = len(E_eps) - 2
chi2_red = chi2 / dof

print(f"chi2_red = {chi2_red:.2f}")

# sortieren für schöneren Plot
ordr = np.argsort(E_eps)
E_eps_s = E_eps[ordr]
eps_s = eps[ordr]
deps_s = deps[ordr]
labels_s = [labels[i] for i in ordr]
residuals_s = residuals[ordr]

E_grid = np.linspace(E_eps.min() * 0.9, E_eps.max() * 1.1, 400)

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    sharex=True,
    figsize=(7, 6),
    gridspec_kw={'height_ratios': [3, 1]}
)

# Daten + Fit
ax1.errorbar(E_eps_s, eps_s, yerr=deps_s, fmt='.', label="Daten")
ax1.plot(E_grid, lin(E_grid, a, b),
         label=f"Fit: ε = aE + b\nχ²_red = {chi2_red:.2f}")

# Beschriftungen
for x, y, label in zip(E_eps_s, eps_s, labels_s):
    ax1.annotate(
        label,
        (x, y),
        textcoords="offset points",
        xytext=(5, 5),
        fontsize=9
    )

ax1.set_ylabel(r"$\varepsilon$")
ax1.set_title("Detektoreffizienz")
ax1.grid(alpha=0.3)
ax1.legend()

# Residuen
ax2.axhline(0, color="black", lw=1)
ax2.errorbar(E_eps_s, residuals_s, yerr=deps_s, fmt='.')
ax2.set_xlabel(r"E$_\gamma$ (keV)")
ax2.set_ylabel("Residuen")
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.show()

'''
# -------------------- Bestimmen der Elektronenmasse --------------------
peak_Na = expected_channel_from_energy(511, 0.476778, 15.559)
x_fit, y_fit, y_err, popt, perr, residuals, chi2_red = fit_peak_with_noise(channels, data_Na, data_NOI, int(peak_Na - 10), int(peak_Na + 10))

A, mu, sigma, B = popt
dA, dmu, dsigma, dB = perr

title = f"Natrium (Fitbereich {int(peak_Na - 10)}-{int(peak_Na + 10)})"
print(f"\n{title}")
print(f"mu    = {mu:.2f} ± {dmu:.2f} (Kanal)")
print(f"sigma = {sigma:.2f} ± {dsigma:.2f} (Kanal)")
print(f"chi2_red = {chi2_red:.2f}")

plot_fit_with_residuals_and_lines(x_fit, y_fit, y_err, popt, residuals, title)

Na_result = {
    "source": source_name,
    "peak_index": i,
    "xmin": xmin, "xmax": xmax,
    "A": A, "dA": dA,
    "mu": mu, "dmu": dmu,
    "sigma": sigma, "dsigma": dsigma,
    "B": B, "dB": dB,
    "chi2_red": chi2_red,
    }

'''