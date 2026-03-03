import numpy as np
from uncertainties import ufloat 
from uncertainties import umath
import pandas as pd
import os
path = os.path.dirname(os.path.abspath(__file__))
os.chdir(path)

a2 = ufloat(1.0, 1.0)
a3 = ufloat(160.0, 1.0)
b1 = ufloat(3.5, 1.0)
b3 = ufloat(3.9, 1.0)
d1 = ufloat(25.0, 1.0)
d2 = ufloat(1.5, 1.0)
    
a1 = np.array([ufloat(19.0, 1.0), ufloat(102.0, 1.0), ufloat(124.5, 1.0), ufloat(131.8, 1.0), ufloat(137.0, 1.0)])
b2 = np.array([ufloat(296.7, 1.0), ufloat(246.7, 1.0), ufloat(226.7, 1.0), ufloat(206.7, 1.0), ufloat(196.7, 1.0)])

a = np.array([a3 - a1[i] - a2 for i in range(len(a1))])
b = np.array([b2[i] - a3 - b1 - b3 for i in range(len(b2))])
d = d1 - d2
theta = np.array([umath.degrees(umath.atan((d/2)/b[i]) + umath.atan((d/2)/a[i])) for i in range(len(a))])

a_std = np.array([a[i].s for i in range(len(a))])
a = np.array([a[i].n for i in range(len(a))])
b_std = np.array([b[i].s for i in range(len(b))])
b = np.array([b[i].n for i in range(len(b))])
theta_std = np.array([theta[i].s for i in range(len(theta))])
theta = np.array([theta[i].n for i in range(len(theta))])

print('d = ',d)

a = a.round(1)
a_std = a_std.round(1)
b = b.round(1)
b_std = b_std.round(1)
theta = theta.round(1)
theta_std = theta_std.round(1)

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
        f"$\\SI{{{a[i]:.1f} \\pm {a_std[i]:.1f}}}{{\\centi\\meter}}$ & "
        f"$\\SI{{{b[i]:.1f} \\pm {b_std[i]:.1f}}}{{\\centi\\meter}}$ & "
        f"$\\SI{{{theta[i]:.1f} \\pm {theta_std[i]:.1f}}}{{\\degree}}$ \\\\"
    )
    latex_rows.append(latex_row)
    latex_rows.append("\\hline")
latex_rows.append("\\end{tabular}")

# Write LaTeX table file
with open('a_b_theta_table.tex', 'w') as f:
    f.write('\n'.join(latex_rows))
print('Generated a_b_theta_table.tex')


