import numpy as np
import matplotlib.pyplot as plt

def plot_fit_with_pull(x_data, data, yerr, func, fit_params, fit_params_unc=None, fitting_range=(0, 1400), xpm = 50, xlabel="channel", ylabel="counts", include_halfwidth_gaussian=False):
    '''
    Wichtig: fitting_range und xpm sind in einheiten der x-Achse zu wählen, nicht zwingend als indizes"
    '''
    x_data = np.array(x_data)
    data = np.array(data)
    yerr = np.array(yerr)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))
    ### Select all x_data points within fitting_range, x_data being array of values can be energy or channel number
    sel_x = (x_data >= fitting_range[0]-xpm) & (x_data <= fitting_range[1]+xpm)

    ax1.errorbar(x_data[sel_x], data[sel_x], yerr=yerr[sel_x], fmt=".", label="Data")
    x_fit = np.linspace(fitting_range[0], fitting_range[1], 1000)
    y_fit = func(x_fit, *fit_params)
    ax1.plot(x_fit, y_fit, label="Fit", color="red")
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    
    if include_halfwidth_gaussian:
        A_fit, mu_fit, sigma_fit = fit_params[0], fit_params[1], fit_params[2]
        halfwidth = sigma_fit * np.sqrt(2 * np.log(2))
        ax1.axvline(mu_fit, color="green", linestyle="--", label="Peak Position")
        ax1.axvline(mu_fit - halfwidth, color="orange", linestyle="--", label="Halfwidth")
        ax1.axvline(mu_fit + halfwidth, color="orange", linestyle="--", label="Halfwidth")
    ax1.legend()
    sel_pull = (x_data >= fitting_range[0]) & (x_data <= fitting_range[1])
    pulls = (data[sel_pull] - func(x_data[sel_pull], *fit_params)) / yerr[sel_pull]
    chi2 = np.sum(pulls**2)
    ndof = len(pulls) - len(fit_params)
    ax2.errorbar(x_data[sel_pull], pulls, yerr=np.ones_like(pulls), fmt=".", label="Pulls")
    ax2.axhline(0, color="red", linestyle="--")
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(r"$\frac{y_i - f(x_i)}{\sigma_i}$")
    plt.tight_layout()
    plt.show()
    return chi2, ndof, chi2/ndof