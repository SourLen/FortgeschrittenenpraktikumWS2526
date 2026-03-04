import numpy as np
from uncertainties import ufloat 
from uncertainties import umath
import pandas as pd
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

E = np.array([ufloat(0.662, 0) / (1 + (0.662/(0.511))*(1-umath.cos(umath.radians(theta[i])))) for i in range(len(theta))])
print('E = ',E)

# calculations for absorbtion

x1 = np.array([umath.sqrt(a[i]**2 + d**2) for i in range(len(a1))])
x2 = np.array([umath.sqrt(b[i]**2 + d**2) for i in range(len(a1))])
xring = umath.sqrt( 2*(d2/2)**2)

density = ufloat(7.8744, 0) # g/cm^3 https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
mac = ufloat(0, 0) # cm^2/g https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z26.html

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