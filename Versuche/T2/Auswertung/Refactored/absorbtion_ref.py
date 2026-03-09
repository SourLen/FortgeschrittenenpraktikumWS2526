import numpy as np
from uncertainties import ufloat, UFloat ### UFloat nur für type hinting in funktionen iwie kp mein vscode meckert bei ufloat
from uncertainties import umath
from scipy.optimize import curve_fit
from ring_geometry_ref import x2, x1, xring
from convential_geometry_ref import *
from scipy.constants import N_A
# calculations for absorbtion

## Data for absorption calculation

### radiation data

E_gamma = ufloat(0.661657, 0.000003) # Mev Anleitung

### material data 

density_al = 2.330E+00 # g/cm^3 https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
E_ref_al = np.array([1.50000E-01, 2.00000E-01, 3.00000E-01, 4.00000E-01, 5.00000E-01, 6.00000E-01, 8.00000E-01]) # MeV https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html
mac_ref_al = np.array([1.378E-01, 1.223E-01, 1.042E-01, 9.276E-02, 8.445E-02, 7.802E-02, 6.841E-02]) # cm^2/g https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html
al_data = [density_al, E_ref_al, mac_ref_al]
atomic_weight_al = 26.981 #g/mol  https://pubchem.ncbi.nlm.nih.gov/#query=Aluminium


density_fe = 7.874E+00 # g/cm^3 https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
E_ref_fe = np.array([1.50000E-01, 2.00000E-01, 3.00000E-01, 4.00000E-01, 5.00000E-01, 6.00000E-01, 8.00000E-01]) # MeV https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z26.html
mac_ref_fe = np.array([1.964E-01, 1.460E-01, 1.099E-01, 9.400E-02, 8.414E-02, 7.704E-02, 6.699E-02]) # cm^2/g https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z26.html
fe_data = [density_fe, E_ref_fe, mac_ref_fe]
atomic_weight_fe = 55.845 #g/mol  https://pubchem.ncbi.nlm.nih.gov/#query=iron
atom_density_al = density_al / atomic_weight_al * N_A # atoms/cm^3
atom_density_fe = density_fe / atomic_weight_fe * N_A # atoms/cm^3

density_air = 1.205E-03 # g/cm^3 https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html
E_ref_air = np.array([1.50000E-01, 2.00000E-01, 3.00000E-01, 4.00000E-01, 5.00000E-01, 6.00000E-01, 8.00000E-01]) # MeV https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
mac_ref_air = np.array([1.356E-01, 1.223E-01, 1.067E-01, 9.549E-02, 8.712E-02, 8.055E-02, 7.074E-02]) # cm^2/g https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
air_data= [density_air, E_ref_air, mac_ref_air]

## functions for absorption calculation

def mac(x, alpha, beta):
    return alpha*x+beta

def get_mac(E, E_ref, mac_ref):
    popt, pcov = curve_fit(mac, E_ref, mac_ref)
    return mac(E, ufloat(popt[0], pcov[0, 0]**0.5), ufloat(popt[1], pcov[1, 1]**0.5))

def get_absorption(E: UFloat, E0: UFloat, x_before:UFloat, x_after: UFloat, x_inside: UFloat, material: tuple, air=air_data)-> UFloat:
    E_steps = np.array([E0, E0, E, E])
    E_ref_steps = np.array([air[1], material[1], material[1], air[1]])
    mac_ref_steps = np.array([air[2], material[2], material[2], air[2]])
    mac_steps = np.array([get_mac(E_steps[i], E_ref_steps[i], mac_ref_steps[i]) for i in range(len(E_steps))])

    density_steps = np.array([air[0], material[0], material[0], air[0]])

    x_steps = np.array((x_before, x_inside, x_inside, x_after))

    eta = 1
    for i in range(len(x_steps)):
        eta *= umath.exp(-mac_steps[i] * density_steps[i] * x_steps[i])
    return eta

E_example = ufloat(0.6, 0.0001)

print('Beispiel: eta = ', get_absorption(E=E_example, E0=E_gamma, x_before=x2[0], x_after=x1[0], x_inside=xring, material=al_data))