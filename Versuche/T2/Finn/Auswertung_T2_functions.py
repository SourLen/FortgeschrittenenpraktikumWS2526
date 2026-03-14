import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# =========================================================
# Grundlegende Modellfunktionen
# =========================================================

def gauss(x, A, mu, sigma, B):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2)) + B

def lin(E, a, b):
    return a * E + b


def inv_lin(ch, m, q):
    return (ch - q) / m


def res_model(E, a, b):
    return np.sqrt(a**2 + (b**2) / E)


# =========================================================
# Datei einlesen
# =========================================================

def read_tka(path, n_channels=1024):
    return np.loadtxt(path)[:n_channels]

# =========================================================
# Physik: Compton-Streuung
# =========================================================

def E_scattered_Cs137(theta_deg):
    E0 = 661.66       # Energie der Cs137 Gamma-Linie
    mec2 = 510.999    # Elektronenruheenergie
    theta = np.deg2rad(theta_deg)
    return E0 / (1 + (E0 / mec2) * (1 - np.cos(theta)))


def expected_channel(theta_deg, m, q): # Erwartete Kanalposition eines Compton-Peaks aus Winkel und Energiekalibration.
    return m * E_scattered_Cs137(theta_deg) + q


def expected_channel_from_energy(E, m, q): # Erwarteter Kanal für eine bekannte Energie.
    return m * E + q


# =========================================================
# Energieumrechnung und Fehler
# =========================================================

def channel_to_energy(ch, m, q): # Rechnet Kanal in Energie um.
    return (ch - q) / m


def channel_to_energy_err(ch, dch, m, dm, q, dq): # Fehlerfortpflanzung für E = (ch - q)/m
    dE_dch = 1 / m
    dE_dq = -1 / m
    dE_dm = -(ch - q) / m**2

    return np.sqrt(
        (dE_dch * dch)**2 +
        (dE_dq * dq)**2 +
        (dE_dm * dm)**2
    )


def theory_energy_err(theta_deg, dtheta_deg): #    Fehler der theoretischen Energie durch Winkelunsicherheit.

    # numerische Ableitung
    h = 1e-3

    dEdtheta = (
        E_scattered_Cs137(theta_deg + h)
        - E_scattered_Cs137(theta_deg - h)
    ) / (2 * h)

    return abs(dEdtheta) * dtheta_deg


# =========================================================
# Gauß-Fits
# =========================================================

def fit_peak_with_noise(x, y, noise, xmin, xmax):
    """
    Gauß-Fit mit vorherigem Abzug einer Rauschmessung.
    """

    mask = (x >= xmin) & (x <= xmax)

    x_fit = x[mask]

    # Untergrund abziehen
    y_fit = y[mask] - noise[mask]

    y_fit[y_fit < 0] = 0

    # Fehler (Poisson Statistik)
    y_err = np.sqrt(y[mask] + noise[mask])
    y_err[y_err == 0] = 1

    # Startwerte
    A0 = np.max(y_fit) - np.min(y_fit)
    mu0 = x_fit[np.argmax(y_fit)]
    sigma0 = 5
    B0 = np.min(y_fit)

    popt, pcov = curve_fit(
        gauss,
        x_fit,
        y_fit,
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


# =========================================================
# Plot: Fit + Residuen
# =========================================================

def plot_fit_with_residuals_and_lines(
    x_fit,
    y_fit,
    y_err,
    popt,
    residuals,
    title,
    x_theo=None,
    theo_label="erwarteter Peak"
):
    """
    Plot des Fits mit Residuen.
    Zusätzlich vertikale Linien für:
    
    - Fit-Peak (rot)
    - theoretische Position (schwarz)
    """

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
        figsize=(6.8, 6.2)
    )

    mu_fit = popt[1]

    ax1.errorbar(x_fit, y_fit, yerr=y_err, fmt=".", label="Daten")
    ax1.plot(x_fit, gauss(x_fit, *popt), label="Gauß-Fit")
    ax1.axvline(
        mu_fit,
        color="red",
        linestyle="--",
        linewidth=1.6,
        label=fr"Fit-Peak $\mu = {mu_fit:.2f}$"
    )

    if x_theo is not None:
        ax1.axvline(
            x_theo,
            color="black",
            linestyle=":",
            linewidth=1.8,
            label=fr"{theo_label} = {x_theo:.2f}"
        )

    ax1.set_ylabel("Counts")
    ax1.set_title(title)
    ax1.grid(alpha=0.3)
    ax1.legend(loc="upper right")
    ax2.axhline(0, lw=1)
    ax2.errorbar(x_fit, residuals, yerr=y_err, fmt=".")
    ax2.set_xlabel("Kanal")
    ax2.set_ylabel("Residuen")
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


# =========================================================
# Hilfsfunktion: Spektren plotten
# =========================================================

def plot_spectra(
    channels,
    spectra,
    noise=None,
    subtract_noise=False,
    figsize=(10, 12),
    sharex=True,
    xlabel="Kanal",
    ylabel=None,
    title_suffix=""
):
    """
    Plottet mehrere Spektren untereinander.

    Parameter
    ---------
    channels : array
        Kanalachse
    spectra : dict
        z.B. {"Cs-137": data_Cs, "Na-22": data_Na, ...}
    noise : array oder None
        Rauschmessung
    subtract_noise : bool
        Falls True, wird noise abgezogen
    """
    names = list(spectra.keys())
    n = len(names)

    fig, axs = plt.subplots(n, 1, figsize=figsize, sharex=sharex)

    if n == 1:
        axs = [axs]

    for ax, name in zip(axs, names):
        y = spectra[name].copy()

        if subtract_noise:
            if noise is None:
                raise ValueError("subtract_noise=True, aber noise=None")
            y = y - noise

        ax.step(channels, y, where="mid", label=name)
        ax.set_title(f"{name}{title_suffix}")
        ax.grid(alpha=0.3)

        if ylabel is not None:
            ax.set_ylabel(ylabel)

        ax.legend(loc="upper right")

    axs[-1].set_xlabel(xlabel)
    plt.tight_layout()
    plt.show()


# =========================================================
# Hilfsfunktion: Spektren mit Peak-Linien plotten
# =========================================================

def plot_spectra_with_peaks(
    channels,
    spectra,
    results,
    peak_energies=None,
    m=None,
    q=None,
    noise=None,
    subtract_noise=False,
    figsize=(10, 12),
    sharex=True,
    xlabel="Kanal",
    ylabel=None
):
    """
    Plottet Spektren und zeichnet Peak-Positionen ein:
    - theoretische Position aus peak_energies + Kalibration
    - gefittete Position mu aus results

    results muss Einträge mit
        "source", "peak_index", "mu"
    enthalten.
    """
    names = list(spectra.keys())
    n = len(names)

    fig, axs = plt.subplots(n, 1, figsize=figsize, sharex=sharex)

    if n == 1:
        axs = [axs]

    for ax, name in zip(axs, names):
        y = spectra[name].copy()

        if subtract_noise:
            if noise is None:
                raise ValueError("subtract_noise=True, aber noise=None")
            y = y - noise

        ax.step(channels, y, where="mid", label="Spektrum")
        ax.set_title(name)
        ax.grid(alpha=0.3)

        if ylabel is not None:
            ax.set_ylabel(ylabel)

        # theoretische Peaks
        theo_drawn = False
        if peak_energies is not None and m is not None and q is not None and name in peak_energies:
            for E in peak_energies[name]:
                mu_theo = expected_channel_from_energy(E, m, q)
                ax.axvline(
                    mu_theo,
                    color="black",
                    linestyle=":",
                    linewidth=1.6,
                    label="theoretischer Peak" if not theo_drawn else None
                )
                theo_drawn = True

        # gefittete Peaks
        fit_drawn = False
        sub_results = [r for r in results if r["source"] == name]
        for r in sub_results:
            ax.axvline(
                r["mu"],
                color="red",
                linestyle="--",
                linewidth=1.6,
                label="Gauß-Fit: μ" if not fit_drawn else None
            )
            fit_drawn = True

        ax.legend(loc="upper right")

    axs[-1].set_xlabel(xlabel)
    plt.tight_layout()
    plt.show()
    
    
# =========================================================
# Automatische Compton-Auswertung
# =========================================================

def fit_all_peaks(files, channels, m, dm, q, dq, dtheta_deg, fit_windows, plot_each_fit=True):

    results = []

    for geometry, geo_data in files.items():
        for angle in sorted(geo_data.keys()):
            entry = geo_data[angle]

            noise = read_tka(entry["NOI"])

            noise_factor = 4 if geometry == "Ring" else 1
            noise_scaled = noise_factor * noise

            for material, path in entry.items():

                if material == "NOI":
                    continue

                y = read_tka(path)

                mu_exp = expected_channel(angle, m, q)

                # Fenster aus Dictionary
                left, right = fit_windows[(geometry, material, angle)]

                xmin = max(0, mu_exp + left)
                xmax = min(len(channels)-1, mu_exp + right)

                x_fit, y_fit, y_err, popt, perr, residuals, chi2_red = fit_peak_with_noise(
                    channels, y, noise_scaled, xmin, xmax
                )

                A, mu, sigma, B = popt
                dA, dmu, dsigma, dB = perr

                E_meas = channel_to_energy(mu, m, q)
                dE_meas = channel_to_energy_err(mu, dmu, m, dm, q, dq)

                E_theo = E_scattered_Cs137(angle)
                dE_theo = theory_energy_err(angle, dtheta_deg)

                results.append({
                    "geometry": geometry,
                    "material": material,
                    "angle": angle,
                    "fit_window": (xmin, xmax),
                    "mu_exp": mu_exp,
                    "mu_fit": mu,
                    "dmu_fit": dmu,
                    "sigma": sigma,
                    "dsigma": dsigma,
                    "A": A,
                    "dA": dA,
                    "E_meas": E_meas,
                    "dE_meas": dE_meas,
                    "E_theo": E_theo,
                    "dE_theo": dE_theo,
                    "chi2_red": chi2_red,
                })

                if plot_each_fit:

                    title = (
                        f"{geometry}, {material}, {angle}°\n"
                        f"$\\mu = {mu:.2f} \\pm {dmu:.2f}$, "
                        f"$E = {E_meas:.2f} \\pm {dE_meas:.2f}$ keV, "
                        f"$\\chi^2_\\mathrm{{red}} = {chi2_red:.2f}$, "
                        f"Fenster: [{xmin:.1f}, {xmax:.1f}]"
                    )

                    plot_fit_with_residuals_and_lines(x_fit, y_fit, y_err, popt, residuals, title, x_theo=mu_exp, theo_label="theor. Peak")
                    
    return results

def plot_energy_vs_angle(results, dtheta_deg):
    fig, ax = plt.subplots(figsize=(8, 5.5))

    combos = sorted(set((r["geometry"], r["material"]) for r in results))

    for geometry, material in combos:
        sub = [r for r in results if r["geometry"] == geometry and r["material"] == material]
        angles = np.array([r["angle"] for r in sub])
        E_meas = np.array([r["E_meas"] for r in sub])
        dE_meas = np.array([r["dE_meas"] for r in sub])

        ax.errorbar(
            angles,
            E_meas,
            yerr=dE_meas,
            fmt="o",
            label=f"{geometry}, {material}"
        )

    angles_theo = np.linspace(
        min(r["angle"] for r in results),
        max(r["angle"] for r in results),
        400
    )
    ax.plot(angles_theo, E_scattered_Cs137(angles_theo), label="Theorie")

    used_angles = sorted(set(r["angle"] for r in results))
    E_theo_pts = np.array([E_scattered_Cs137(a) for a in used_angles])
    dE_theo_pts = np.array([theory_energy_err(a, dtheta_deg) for a in used_angles])

    ax.errorbar(
        used_angles,
        E_theo_pts,
        yerr=dE_theo_pts,
        fmt="s",
        capsize=3,
        label="Theorie an Messwinkeln"
    )

    ax.set_xlabel("Winkel / °")
    ax.set_ylabel("Peak-Energie / keV")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.show()
    
def compare_geometries_same_angle(results, angle=50):
    sub = [r for r in results if r["angle"] == angle]

    print(f"\nVergleich bei {angle}°:")
    for r in sub:
        print(
            f'{r["geometry"]:13s} {r["material"]:2s}: '
            f'E = {r["E_meas"]:.2f} ± {r["dE_meas"]:.2f} keV'
        )

    # Ring-Al mit Konventionell-Al vergleichen
    ring_al = next(r for r in sub if r["geometry"] == "Ring" and r["material"] == "Al")
    conv_al = next(r for r in sub if r["geometry"] == "Konventionell" and r["material"] == "Al")

    dE = ring_al["E_meas"] - conv_al["E_meas"]
    dE_err = np.sqrt(ring_al["dE_meas"]**2 + conv_al["dE_meas"]**2)
    nsigma = dE / dE_err

    print("\nRing-Al vs. Konventionell-Al:")
    print(f"ΔE = {dE:.2f} ± {dE_err:.2f} keV  ({nsigma:.2f}σ)")

    # optional: Fe dazu
    conv_fe = next(r for r in sub if r["geometry"] == "Konventionell" and r["material"] == "Fe")
    dE_fe = ring_al["E_meas"] - conv_fe["E_meas"]
    dE_fe_err = np.sqrt(ring_al["dE_meas"]**2 + conv_fe["dE_meas"]**2)
    nsigma_fe = dE_fe / dE_fe_err

    print("Ring-Al vs. Konventionell-Fe:")
    print(f"ΔE = {dE_fe:.2f} ± {dE_fe_err:.2f} keV  ({nsigma_fe:.2f}σ)")