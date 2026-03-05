import numpy as np
from uncertainties import ufloat 
from uncertainties import umath
import pandas as pd
from scipy.optimize import curve_fit
import os
path = os.path.dirname(os.path.abspath(__file__))
os.chdir(path)

#measurements

a2 = ufloat(1.0, 0.1/np.sqrt(12))
a3 = ufloat(160.0, 0.3)
b1 = ufloat(3.5, 0.1/np.sqrt(12))
b3 = ufloat(3.9, 0.1/np.sqrt(12))
d1 = ufloat(25.0, 0.1/np.sqrt(12))
d2 = ufloat(1.5, 0.1/np.sqrt(12))
    
a1 = np.array([ufloat(19.0, 0.3), ufloat(102.0, 0.3), ufloat(124.5, 0.3), ufloat(131.8, 0.3), ufloat(137.0, 0.3)])
b2 = np.array([ufloat(296.7, 0.3), ufloat(246.7, 0.3), ufloat(226.7, 0.3), ufloat(206.7, 0.3), ufloat(196.7, 0.3)])

# calculations for a, b and theta

a = np.array([a3 - a1[i] - a2 for i in range(len(a1))])
b = np.array([b2[i] - a3 - b1 - b3 for i in range(len(b2))])
d = d1 - d2
theta = np.array([umath.degrees(umath.atan((d/2)/b[i]) + umath.atan((d/2)/a[i])) for i in range(len(a))])

print('d = ',d)

# Energien der Compton gestreuten Photonen

E_diffracted = np.array([ufloat(0.662, 0) / (1 + (0.662/(0.511))*(1-umath.cos(umath.radians(theta[i])))) for i in range(len(theta))])

# calculations for absorbtion

## Data for absorption calculation

x1 = np.array([umath.sqrt(a[i]**2 + (d/2)**2) for i in range(len(a1))])
x2 = np.array([umath.sqrt(b[i]**2 + (d/2)**2) for i in range(len(a1))])
xring = d2*umath.sqrt(1/2)

density_al = 2.330E+00 # g/cm^3 https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
E_ref_al = np.array([1.50000E-01, 2.00000E-01, 3.00000E-01, 4.00000E-01, 5.00000E-01, 6.00000E-01, 8.00000E-01]) # MeV https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html
mac_ref_al = np.array([1.378E-01, 1.223E-01, 1.042E-01, 9.276E-02, 8.445E-02, 7.802E-02, 6.841E-02]) # cm^2/g https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html
density_air = 1.205E-03 # g/cm^3 https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html
E_ref_air = np.array([1.50000E-01, 2.00000E-01, 3.00000E-01, 4.00000E-01, 5.00000E-01, 6.00000E-01, 8.00000E-01]) # MeV https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
mac_ref_air = np.array([1.356E-01, 1.223E-01, 1.067E-01, 9.549E-02, 8.712E-02, 8.055E-02, 7.074E-02]) # cm^2/g https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html


## functions for absorption calculation

def mac(E, alpha, beta):
    return alpha*E+beta

def get_mac(E, E_ref, mac_ref):
    popt, pcov = curve_fit(mac, E_ref, mac_ref)
    return mac(E, ufloat(popt[0], pcov[0, 0]**0.5), ufloat(popt[1], pcov[1, 1]**0.5))

def get_absorption(E, E0, x_before, x_after, x_inside, material):
    x_steps = np.array((x_before, x_inside, x_inside, x_after))
    E_steps = np.array([E0, E0, E, E])
    density_steps = np.array([density_air, material[0], material[0], density_air])
    E_ref_steps = np.array([E_ref_air, material[1], material[1], E_ref_air])
    mac_ref_steps = np.array([mac_ref_air, material[2], material[2], mac_ref_air])
    mac_steps = np.array([get_mac(E_steps[i], E_ref_steps[i], mac_ref_steps[i]) for i in range(len(E_steps))])
    eta = 1
    for i in range(len(x_steps)):
        eta *= umath.exp(-mac_steps[i] * density_steps[i] * x_steps[i])
    return eta

al = [density_al, E_ref_al, mac_ref_al]

print(get_absorption(E_diffracted[0], 0.662, x2[0], x1[0], xring, al))


# Generate LaTeX table rows
latex_rows = []
latex_rows.append("\\begin{tabular}{|c|c|c|}")
latex_rows.append("\\hline")
latex_row = (
    f"$a$ & "
    f"$b$ & "
    f"$\\theta_{{\\text{{Messung}}}}$ \\\\"
)
latex_rows.append(latex_row)
latex_rows.append("\\hline")
for i in range(len(a)):
    latex_row = (
        f"$\\SI{{{a[i].n:.2f} \\pm {a[i].s:.2f}}}{{\\centi\\meter}}$ & "
        f"$\\SI{{{b[i].n:.2f} \\pm {b[i].s:.2f}}}{{\\centi\\meter}}$ & "
        f"$\\SI{{{theta[i].n:.2f} \\pm {theta[i].s:.2f}}}{{\\degree}}$ \\\\"
    )
    latex_rows.append(latex_row)
    latex_rows.append("\\hline")
latex_rows.append("\\end{tabular}")

# Write LaTeX table file
with open('a_b_theta_table.tex', 'w') as f:
    f.write('\n'.join(latex_rows))
print('Generated a_b_theta_table.tex')