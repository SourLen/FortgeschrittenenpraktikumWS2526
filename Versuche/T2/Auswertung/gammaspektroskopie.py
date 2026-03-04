import scipy.optimize as opt
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from plotting_functions import plot_fit_with_pull
from activities import SOURCES
# Einlesen der Daten


### Todo: Eine Funktion fürs einlesen, damit nicht der ganze code 5mal steht

def load_data(filename, noise=None, alpha_noise=None):
    data = np.loadtxt(filename)[2:]
    if noise is not None and alpha_noise is not None:   
        data = data - alpha_noise * noise  # Skip the first two lines
    return data
Noise = load_data("../Messdaten/Gammaspektroskopie/noise.tka", None, None)
Cs137 = load_data("../Messdaten/Gammaspektroskopie/Cs137_MH851_3min.tka", Noise, 3/5)
Co60 = load_data("../Messdaten/Gammaspektroskopie/Co60_LP213_3min.tka", Noise, 3/5)
Eu152 = load_data("../Messdaten/Gammaspektroskopie/Eu152_MH850_3min.tka", Noise, 3/5)
Na22 = load_data("../Messdaten/Gammaspektroskopie/Na22_MH852_3min.tka", Noise, 3/5)

def calculate_unc_after_noise
cs_uncertainty = np.sqrt(Cs137 + (3/5)**2 * Noise)
co_uncertainty = np.sqrt(Co60 + (3/5)**2 * Noise)
eu_uncertainty = np.sqrt(Eu152 + (3/5)**2 * Noise)
na_uncertainty = np.sqrt(Na22 + (3/5)**2 * Noise)
fig, ax = plt.subplots(4, 1, figsize=(10, 15), sharex=False)
ax[0].scatter(np.arange(len(cs_no_noise)), cs_no_noise, label="Cs-137", marker=".", s=6)
ax[1].scatter(np.arange(len(co_no_noise)), co_no_noise, label="Co-60", marker = ".", s=6)
ax[2].scatter(np.arange(len(eu_no_noise)), eu_no_noise, label="Eu-152", marker = ".", s=6)
ax[3].scatter(np.arange(len(na_no_noise)), na_no_noise, label="Na-22", marker = ".", s=6)
### Add x and y axis labels + legend
for axis in ax:
    axis.set_ylabel("counts")
    axis.legend()
plt.xlabel("channel")
plt.tight_layout()
plt.show()
plt.close()



### Peak-Fitting

def find_channel_fit(data, unc, fitting_range, plot=False):
    def gaussian_withbg(x, a, x0, sigma, p_0, p_1):
        return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + (p_0 + p_1 * x)
    x_data = np.arange(len(data))
    unc = np.where(unc == 0, 1, unc)  # Avoid 0 uncertainties
    popt, pcov = curve_fit(gaussian_withbg, x_data[fitting_range[0]:fitting_range[1]], data[fitting_range[0]:fitting_range[1]],
                           p0=[max(data[fitting_range[0]:fitting_range[1]]), np.mean(fitting_range), 10, 0, 0], 
                           sigma=unc[fitting_range[0]:fitting_range[1]], absolute_sigma=True, maxfev = 20000)
    
    A_fit, mu_fit, sigma_fit, p_0_fit, p_1_fit = popt
    perr = np.sqrt(np.diag(pcov))
    if plot:
        plot_fit_with_pull(x_data, data, unc, gaussian_withbg, popt, perr, fitting_range, include_halfwidth_gaussian=True)
        
    N_peak = ufloat(A_fit, perr[0]) * ufloat(sigma_fit, perr[2]) * np.sqrt(2 * np.pi)  # type: ignore ### Area under the Gaussian, proportional to the number of counts in the peak
    return mu_fit, perr[1], ufloat(sigma_fit, perr[2]), N_peak

fitting_ranges = [(450, 550), (835, 905), (950, 1020), (90, 120), (170 , 215),(240, 290), (545, 620), (755, 860)]
peak_1_cs, peak_1_cs_unc, u_sig_1_cs, N_peak_1_cs = find_channel_fit(cs_no_noise, cs_uncertainty, fitting_ranges[0])
print(f"Cs-137 peak: {peak_1_cs:.2f} ± {peak_1_cs_unc:.2f} channels")

peak1_co, peak_1_co_unc, usig_1_co, N_peak_1_co = find_channel_fit(co_no_noise, co_uncertainty, fitting_ranges[1])
peak2_co, peak_2_co_unc, usig_2_co, N_peak_2_co = find_channel_fit(co_no_noise, co_uncertainty, fitting_ranges[2])
print(f"Co-60 peaks: {peak1_co:.2f} ± {peak_1_co_unc:.2f} channels, {peak2_co:.2f} ± {peak_2_co_unc:.2f} channels")
### Todo: siehe oben
peak_1_eu, peak_1_eu_unc, u_sig_1_eu, N_peak_1_eu = find_channel_fit(eu_no_noise, eu_uncertainty, fitting_ranges[3])
peak_2_eu, peak_2_eu_unc, u_sig_2_eu, N_peak_2_eu = find_channel_fit(eu_no_noise, eu_uncertainty, fitting_ranges[4])
peak_3_eu, peak_3_eu_unc, u_sig_3_eu, N_peak_3_eu = find_channel_fit(eu_no_noise, eu_uncertainty, fitting_ranges[5])
peak_4_eu, peak_4_eu_unc, u_sig_4_eu, N_peak_4_eu = find_channel_fit(eu_no_noise, eu_uncertainty, fitting_ranges[6])
peak_5_eu, peak_5_eu_unc, u_sig_5_eu, N_peak_5_eu = find_channel_fit(eu_no_noise, eu_uncertainty, fitting_ranges[7], plot=True)
print(f"Eu-152 peaks: {peak_1_eu:.2f} ± {peak_1_eu_unc:.2f} channels, {peak_2_eu:.2f} ± {peak_2_eu_unc:.2f} channels, {peak_3_eu:.2f} ± {peak_3_eu_unc:.2f} channels, {peak_4_eu:.2f} ± {peak_4_eu_unc:.2f} channels")

energies = np.array([661.66, 1173.23, 1332.49, 121.78, 244.7, 344.28, 778.9, 1085.9])
channels = np.array([peak_1_cs, peak1_co, peak2_co, peak_1_eu, peak_2_eu, peak_3_eu, peak_4_eu, peak_5_eu])
u_channels = np.array([peak_1_cs_unc, peak_1_co_unc, peak_2_co_unc, peak_1_eu_unc, peak_2_eu_unc, peak_3_eu_unc, peak_4_eu_unc, peak_5_eu_unc])
def linear(x, m, b):
    return m * x + b
popt, pcov = curve_fit(linear, energies, channels, sigma=np.array([peak_1_co_unc, peak_1_co_unc, peak_2_co_unc, peak_1_eu_unc, peak_2_eu_unc, peak_3_eu_unc, peak_4_eu_unc, peak_5_eu_unc]), absolute_sigma=True)
m_fit, b_fit = popt
perr = np.sqrt(np.diag(pcov))
_,_,chi2ndf = plot_fit_with_pull(energies, channels, u_channels, linear, popt, perr, fitting_range=(0, 1400), xlabel="Energy (keV)", ylabel="Channel")
print("Chi2/ndf for energy calibration fit:", chi2ndf)
### Energieauflösung
sigma_channels = [u_sig_1_cs, usig_1_co, usig_2_co, u_sig_1_eu, u_sig_2_eu, u_sig_3_eu, u_sig_4_eu, u_sig_5_eu]
u_m_fit = ufloat(m_fit, perr[0]) ### Aus der linearen Beziehung Channel = m*E_gamma + b
u_b_fit = ufloat(b_fit, perr[1])
u_delta_energies = [2*np.sqrt(np.log(4))*sigma_channels[i] / u_m_fit for i in range(len(sigma_channels))] ### FWHM = 2*sqrt(ln(4))*sigma, umgerechnet in Energie durch Division durch m
u_energies_peak_fit = []
for i in range(len(channels)):
    u_channel = ufloat(channels[i], u_channels[i])
    energy = (u_channel - u_b_fit) / u_m_fit # type: ignore #### Korrelation? Gottweiß, evtl per Hand
    u_energies_peak_fit.append(energy)

### Plotting (delta E/ E)^2 vs 1/E
delta_E_over_E_squared = [(u_delta_energies[i] / u_energies_peak_fit[i])**2 for i in range(len(u_delta_energies))]
inverse_energies = [1/u_energies_peak_fit[i] for i in range(len(u_energies_peak_fit))]

plt.errorbar([i.n for i in inverse_energies], [i.n for i in delta_E_over_E_squared], yerr=[i.s for i in delta_E_over_E_squared], fmt="o")
plt.xlabel("1/E (keV⁻¹)")
plt.ylabel("(ΔE/E)²")
plt.title("Energy Resolution")
plt.show()

#### Alternative: Direct fit onto Delta E/E vs E, with f(E) = sqrt(a^2+b^2/E)
def fit_func(E, a, b):
    return np.sqrt(a**2 + b**2 / E)
y_data = [u_delta_energies[i] / u_energies_peak_fit[i] for i in range(len(u_delta_energies))]
for i in range(len(y_data)):
    print(f"deltaE/E = {y_data[i]:.4f} for E = {u_energies_peak_fit[i]:.2f} keV")
x_data = [u_energies_peak_fit[i].n for i in range(len(u_energies_peak_fit))]


popt, pcov = curve_fit(fit_func, x_data, [y.n for y in y_data], sigma=[y.s for y in y_data], absolute_sigma=True)
a_fit, b_fit = popt
perr = np.sqrt(np.diag(pcov))
u_a_fit, u_b_fit = ufloat(a_fit, perr[0]), ufloat(b_fit, perr[1])
print(u_a_fit, u_b_fit)
_, _, chi2ndof =plot_fit_with_pull(x_data, [y.n for y in y_data], [y.s for y in y_data], fit_func, popt, perr, fitting_range=(100, 1400), xlabel="Energy (keV)", ylabel=r"$\Delta E / E$")
### Scheiße
print(chi2ndof)
### Efficiency Calculation.
N_peaks = [N_peak_1_cs, N_peak_1_co, N_peak_2_co, N_peak_1_eu, N_peak_2_eu, N_peak_3_eu, N_peak_4_eu, N_peak_5_eu]
Time = 180  # seconds
activity_dict = {source.name: source.activity_on() for source in SOURCES.values()}
I_gamma = []
