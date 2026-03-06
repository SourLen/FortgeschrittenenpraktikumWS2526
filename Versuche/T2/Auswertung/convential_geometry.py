import numpy as np
from uncertainties import ufloat 
from uncertainties import umath
import pandas as pd

# measurements

rT = ufloat(6.7, 0.01/np.sqrt(12))
s0 = ufloat(0.78, 0.01/np.sqrt(12))
s1 = ufloat(0.675, 0.01/np.sqrt(12))
s2 = ufloat(5.075, 0.01/np.sqrt(12))
s3 = ufloat(9.350, 0.01/np.sqrt(12))
Ds = ufloat(1.1, 0.01/np.sqrt(12))
Dk = ufloat(2.445, 0.01/np.sqrt(12))
Dz = ufloat(2.445, 0.01/np.sqrt(12))

# distances

x_luft1 = rT - s0
x_luft2 = rT + s1 + s2 + s3
x_material = Ds/2