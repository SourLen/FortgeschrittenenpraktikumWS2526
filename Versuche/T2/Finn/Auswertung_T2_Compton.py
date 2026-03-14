import numpy as np
import matplotlib.pyplot as plt
from Auswertung_T2_functions import *

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "figure.titlesize": 16,
})

channels = np.arange(1024)

m = 0.476778
dm = 0.000444
q = 15.559
dq = 0.166
dtheta_deg = 0.5 / np.sqrt(12)

fit_windows = {

    ("Ring","Al",10): (-25,25),
    ("Ring","Al",20): (-30,40),
    ("Ring","Al",30): (-40,50),
    ("Ring","Al",40): (-40,50),
    ("Ring","Al",50): (-40,50),

    ("Konventionell","Al",50): (-60,40),
    ("Konventionell","Fe",50): (-60,40),

    ("Konventionell","Al",80): (-30,40),
    ("Konventionell","Fe",80): (-30,40),

    ("Konventionell","Al",90): (-30,30),
    ("Konventionell","Fe",90): (-30,30),

    ("Konventionell","Al",105): (-20,30),
    ("Konventionell","Fe",105): (-20,30),

    ("Konventionell","Al",135): (-20,30),
    ("Konventionell","Fe",135): (-20,30),
}

files = {
    "Ring": {
        10: {
            "Al": "01_T2_Daten/02_Compton/Ring/01_10deg_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Ring/02_10deg_rausch_5min.tka",
        },
        20: {
            "Al": "01_T2_Daten/02_Compton/Ring/03_20deg_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Ring/04_20deg_rausch_5min.tka",
        },
        30: {
            "Al": "01_T2_Daten/02_Compton/Ring/05_30deg_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Ring/06_30deg_rausch_5min.tka",
        },
        40: {
            "Al": "01_T2_Daten/02_Compton/Ring/07_40deg_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Ring/08_40deg_rausch_5min.tka",
        },
        50: {
            "Al": "01_T2_Daten/02_Compton/Ring/09_50deg_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Ring/10_50deg_rausch_5min.tka",
        },
    },
    "Konventionell": {
        50: {
            "Al": "01_T2_Daten/02_Compton/Konventionell/02_50deg_Al_20min.tka",
            "Fe": "01_T2_Daten/02_Compton/Konventionell/01_50deg_Fe_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Konventionell/00_50deg_rausch_20min.tka",
        },
        80: {
            "Al": "01_T2_Daten/02_Compton/Konventionell/04_80deg_Al_20min.tka",
            "Fe": "01_T2_Daten/02_Compton/Konventionell/03_80deg_Fe_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Konventionell/05_80deg_rausch_20min.tka",
        },
        90: {
            "Al": "01_T2_Daten/02_Compton/Konventionell/07_90deg_Al_20min.tka",
            "Fe": "01_T2_Daten/02_Compton/Konventionell/06_90deg_Fe_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Konventionell/08_90deg_rausch_20min.tka",
        },
        105: {
            "Al": "01_T2_Daten/02_Compton/Konventionell/10_105deg_Al_20min.tka",
            "Fe": "01_T2_Daten/02_Compton/Konventionell/09_105deg_Fe_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Konventionell/11_105deg_rausch_20min.tka",
        },
        135: {
            "Al": "01_T2_Daten/02_Compton/Konventionell/13_135deg_Al_20min.tka",
            "Fe": "01_T2_Daten/02_Compton/Konventionell/12_135deg_Fe_20min.tka",
            "NOI": "01_T2_Daten/02_Compton/Konventionell/14_135deg_rausch_20min.tka",
        },
    }
}



# ------------- Aufgabe 2 -------------
results = fit_all_peaks(
    files=files,
    channels=channels,
    m = m, dm = dm, q = q, dq = dq,
    dtheta_deg = dtheta_deg,
    fit_windows = fit_windows,
    plot_each_fit=False)

#plot_energy_vs_angle(results, dtheta_deg)

# Ausgeben aller wichtigen Daten:
for r in results:
    xmin, xmax = r["fit_window"]
    #print(
     ##   f'{r["geometry"]:13s}  {r["material"]:2s}  {r["angle"]:>3}°   '
       # f'Fenster = ({xmin:.1f}, {xmax:.1f})   '
        #f'mu = {r["mu_fit"]:7.2f} ± {r["dmu_fit"]:5.2f}   '
        #f'E = {r["E_meas"]:7.2f} ± {r["dE_meas"]:5.2f} keV   '
        #f'chi2_red = {r["chi2_red"]:.2f}'
    #)
    
# ------------- Aufgabe 1 -------------
compare_geometries_same_angle(results, angle=50)

angle = 50

ring_al = read_tka(files["Ring"][angle]["Al"])
ring_noi = 4 * read_tka(files["Ring"][angle]["NOI"])

conv_al = read_tka(files["Konventionell"][angle]["Al"])
conv_fe = read_tka(files["Konventionell"][angle]["Fe"])
conv_noi = read_tka(files["Konventionell"][angle]["NOI"])

mu_theo = expected_channel(angle, m, q)

sub = [r for r in results if r["angle"] == 50]

mu_ring_al = next(r["mu_fit"] for r in sub if r["geometry"] == "Ring" and r["material"] == "Al")
mu_conv_al = next(r["mu_fit"] for r in sub if r["geometry"] == "Konventionell" and r["material"] == "Al")
mu_conv_fe = next(r["mu_fit"] for r in sub if r["geometry"] == "Konventionell" and r["material"] == "Fe")

ch = channels[:512]

fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

axs[0].plot(ch, (ring_al - ring_noi)[:512], lw=1.2, label="Ring, Al")
axs[1].plot(ch, (conv_al - conv_noi)[:512], lw=1.2, label="Konventionell, Al")
axs[2].plot(ch, (conv_fe - conv_noi)[:512], lw=1.2, label="Konventionell, Fe")

# Vertikale Linien auf ALLE Achsen zeichnen
for ax in axs:
    ax.axvline(mu_theo, color="black", linestyle=":", linewidth=1.8,
               label=f"theor. Peak = {mu_theo:.1f}")
    ax.axvline(mu_ring_al, color="red", linestyle="--", linewidth=1.5,
               label=f"Fit Ring-Al = {mu_ring_al:.1f}")
    ax.axvline(mu_conv_al, color="blue", linestyle="--", linewidth=1.5,
               label=f"Fit Konv.-Al = {mu_conv_al:.1f}")
    ax.axvline(mu_conv_fe, color="green", linestyle="--", linewidth=1.5,
               label=f"Fit Konv.-Fe = {mu_conv_fe:.1f}")

    ax.set_ylabel("Counts")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 512)

axs[0].set_title("Ring, Al, 50°")
axs[1].set_title("Konventionell, Al, 50°")
axs[2].set_title("Konventionell, Fe, 50°")
axs[2].set_xlabel("Kanal")

# Legende nur oben einmal
axs[0].legend(loc="upper right")

plt.tight_layout()
#plt.show()

# ------------- Aufgabe 3 -------------
# Geometrische Daten:
# für 50 grad wurde Ring verwendet
r0_angle = {10: 1.30, 20: 0.79, 30: 0.60, 40: 0.39, 50: 0.29, 80: 0.0495, 90: 0.0495, 105: 0.0495, 135: 0.0495}  # m
dr0_angle = {10: 0.005, 20: 0.005, 30: 0.005, 40: 0.005, 50: 0.005, 80: 0.0005, 90: 0.0005, 105: 0.0005, 135: 0.0005}  # m
r_angle = {10: 1.40, 20: 0.57, 30: 0.34, 40: 0.27, 50: 0.22, 80: 0.1882, 90: 0.1882, 105: 0.1882, 135: 0.1882}  # m
dr_angle = {10: 0.005, 20: 0.005, 30: 0.005, 40: 0.005, 50: 0.005, 80: 0.0015, 90: 0.0015, 105: 0.0015, 135: 0.0015}  # m
F_d = {10: 5.13e-3, 20: 5.13e-3, 30: 5.13e-3, 40: 5.13e-3, 50: 5.13e-3, 80: 1.26e-3, 90: 1.26e-3, 105: 1.26e-3, 135: 1.26e-3}  # m^2 --- 4.04**2*np.pi bzw. 2**2*np.pi -> noch in cm**2
dF_d = {10: 0.00013, 20: 0.00013, 30: 0.00013, 40: 0.00013, 50: 0.00013, 80: 0.00006, 90: 0.00006, 105: 0.00006, 135: 0.00006}  # m^2 --- 2*np.pi*4.04*0.05 bzw.2*np.pi*2*0.05 -> noch in cm**2

# andere Daten
# Anzahl der Elektronen im Zylinder bzw. Ring
N_Fe_Kon = 5.15 * 1e24
dN_Fe_Kon = 0.04 * 1e24
N_Al_Kon = 1.891 * 1e24
dN_Al_Kon =  0.0006 * 1e24
N_Al_Ring = 6.17 * 1e25
dN_Al_Ring =  0.0006 * 1e24

# --------------------------------------------------
# Elektronenzahl pro Winkel
# hier nur Aluminium
# --------------------------------------------------
N_angle = {10: N_Al_Ring, 20: N_Al_Ring, 30: N_Al_Ring, 40: N_Al_Ring, 50: N_Al_Ring, 80: N_Al_Kon, 90: N_Al_Kon, 105: N_Al_Kon, 135: N_Al_Kon}
dN_angle = {10: dN_Al_Ring, 20: dN_Al_Ring, 30: dN_Al_Ring, 40: dN_Al_Ring, 50: dN_Al_Ring, 80:dN_Al_Kon, 90: dN_Al_Kon, 105: dN_Al_Kon, 135: dN_Al_Kon}

I = 0.85
dI = 0.002  

activity = 1.660 * 10**7 # in Bq
dactivity = 5 * 10**4

# Bestimmen der Zählraten bei den Photopeaks
t = 1200


M_angle = {}
dM_angle = {}
'''


    if angle == 50:   # gewichtetes Mittel aus Ring + Konventionell
        sigma = np.array([v["sigma"] for v in vals])
        dsigma = np.array([v["dsigma"] for v in vals])
        
        w = 1 / dsigma**2
        sigma_mean = np.sum(w * sigma) / np.sum(w)
        dsigma_mean = np.sqrt(1 / np.sum(w))

        A = np.array([v["A"] for v in vals])
        dA = np.array([v["dA"] for v in vals])
        
        w_A = 1 / dA**2
        A_mean = np.sum(w_A * A) / np.sum(w_A)
        dA_mean = np.sqrt(1 / np.sum(w_A))
        
'''
for angle in sorted(set(r["angle"] for r in results)):

    vals = [r for r in results if r["angle"] == angle and r["material"] == "Al"]


    if 1 < 0:
        pass
    else:   # nur eine Messung
        sigma_mean = vals[0]["sigma"]
        dsigma_mean = vals[0]["dsigma"]
        
        A_mean = vals[0]["A"]
        dA_mean = vals[0]["dA"]
        
    M =  np.sqrt(np.pi * 2) *  sigma_mean * A_mean / t
    dM = M * np.sqrt((dA_mean/A_mean)**2 + (dsigma_mean / sigma_mean)**2)
    M_angle[angle] = M
    dM_angle[angle] = dM

# --------------------------------------------------
# Platzhalter für Detektoreffizienz epsilon(E)
# -> später einfach durch deine gefittete Funktion ersetzen
# --------------------------------------------------
def epsilon_det(E):
    return -1/2500 * E + 0.754

def depsilon_det(E):
    # Beispiel / Platzhalter
    return 0.0

# --------------------------------------------------
# Platzhalter für Schwächungskoeffizienten mu(E)
# -> später mit echten Werten/Fits ersetzen
# --------------------------------------------------
def mu_air(E):
    return 0.0

def dmu_air(E):
    return 0.0

def mu_al(E):
    return 0.0

def dmu_al(E):
    return 0.0

# --------------------------------------------------
# Wegstrecken für Absorption
# Hier musst du später deine abgeschätzten/geometrischen Werte einsetzen
# in Metern
# --------------------------------------------------
R_Al = 0.00600       # m
dR_Al = 0.000015     # m

x_air_before = {angle: r0_angle[angle] - R_Al for angle in r0_angle}
dx_air_before = {angle: np.sqrt(dr0_angle[angle]**2 + dR_Al**2) for angle in r0_angle}

x_al_before = {angle: R_Al for angle in r0_angle}
dx_al_before = {angle: dR_Al for angle in r0_angle}

x_air_after = {angle: r_angle[angle] - R_Al for angle in r_angle}
dx_air_after = {angle: np.sqrt(dr_angle[angle]**2 + dR_Al**2) for angle in r_angle}

x_al_after = {angle: R_Al for angle in r_angle}
dx_al_after = {angle: dR_Al for angle in r_angle}

# --------------------------------------------------
# Absorptionsfaktor eta
# --------------------------------------------------
E_in = 661.66 
eta_angle = {}
deta_angle = {}

for angle in sorted(M_angle.keys()):
    E_out = E_scattered_Cs137(angle)

    eta = (
        np.exp(-mu_air(E_in)  * x_air_before[angle]) *
        np.exp(-mu_al(E_in)   * x_al_before[angle])  *
        np.exp(-mu_al(E_out)  * x_al_after[angle])   *
        np.exp(-mu_air(E_out) * x_air_after[angle])
    )

    # logarithmische Fehlerfortpflanzung
    rel2 = (
        (x_air_before[angle] * dmu_air(E_in))**2 +
        (mu_air(E_in) * dx_air_before[angle])**2 +
        (x_al_before[angle] * dmu_al(E_in))**2 +
        (mu_al(E_in) * dx_al_before[angle])**2 +
        (x_al_after[angle] * dmu_al(E_out))**2 +
        (mu_al(E_out) * dx_al_after[angle])**2 +
        (x_air_after[angle] * dmu_air(E_out))**2 +
        (mu_air(E_out) * dx_air_after[angle])**2
    )

    deta = eta * np.sqrt(rel2)

    eta_angle[angle] = eta
    deta_angle[angle] = deta

    

# --------------------------------------------------
# differentieller Wirkungsquerschnitt
# --------------------------------------------------
dsdo_angle = {}
ddsdo_angle = {}
ind = 0
import pandas as pd
data = {"angle": [],
        "material": [],
        "geometry": [],
        "r_0^2": [],
        "r^2": [],
        "rate M": [],
        "activity": [],
        "I": [],
        "efficiency": [],
        "eta": [],
        "n_electrons": [],
        "F_D_cm2": [],
        "dsdo_cm2": []
        }
for angle in sorted(M_angle.keys()):
    E_out = E_scattered_Cs137(angle)

    eps = epsilon_det(E_out)
    deps = depsilon_det(E_out)

    eta = eta_angle[angle]
    deta = deta_angle[angle]
    M = M_angle[angle]
    dM = dM_angle[angle]

    r0 = r0_angle[angle]
    dr0 = dr0_angle[angle]

    r = r_angle[angle]
    dr = dr_angle[angle]

    FD = F_d[angle]
    dFD = dF_d[angle]

    N = N_angle[angle]
    dN = dN_angle[angle]

    dsdo = M * 4 * np.pi * r0**2 * r**2 / (activity * I * eta * eps * N * FD)

    rel2 = (
        (dM / M)**2 +
        (2 * dr0 / r0)**2 +
        (2 * dr / r)**2 +
        (dactivity / activity)**2 +
        (dI / I)**2 +
        (deta / eta)**2 +
        (deps / eps)**2 +
        (dN / N)**2 +
        (dFD / FD)**2
    )

    ddsdo = dsdo * np.sqrt(rel2)

    dsdo_angle[angle] = dsdo
    ddsdo_angle[angle] = ddsdo
    
    data["angle"].append(angle)
    data["material"].append("Al")
    data["geometry"].append("-")
    data["r_0^2"].append(r0**2*1e4)
    data["r^2"].append(r**2*1e4)
    data["rate M"].append(M)
    data["activity"].append(activity)
    data["I"].append(I)
    data["efficiency"].append(eps)
    data["eta"].append(eta)
    data["n_electrons"].append(N)
    data["F_D_cm2"].append(FD*1e4)
    data["dsdo_cm2"].append(dsdo*1e4)
    ind += 1
    print(FD/r**2)
pd.DataFrame(data).to_csv("compton_scattering_results_finn.csv", sep=";", decimal=",", index=False, encoding="utf-8")
angles = np.array(sorted(dsdo_angle.keys()))
dsdo_vals = np.array([dsdo_angle[a] for a in angles])
ddsdo_vals = np.array([ddsdo_angle[a] for a in angles])

plt.figure()
plt.errorbar(angles, dsdo_vals, yerr=ddsdo_vals, fmt='.')
plt.xlabel("Streuwinkel (°)")
plt.ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\Omega$ (m$^2$/sr)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

def klein_nishina(theta_deg):
    r_e = 2.8179403262e-15  # m
    theta = np.deg2rad(theta_deg)

    E_ratio = E_scattered_Cs137(theta_deg) / E_in

    return 0.5 * r_e**2 * E_ratio**2 * (E_ratio + 1/E_ratio - np.sin(theta)**2)

theta_grid = np.linspace(min(angles), max(angles), 400)

plt.figure()
plt.errorbar(angles, dsdo_vals, yerr=ddsdo_vals, fmt='.', label="Messung")
plt.plot(theta_grid, klein_nishina(theta_grid), label="Klein-Nishina")
plt.xlabel("Streuwinkel (°)")
plt.ylabel(r"$\mathrm{d}\sigma/\mathrm{d}\Omega$ (m$^2$/sr)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
