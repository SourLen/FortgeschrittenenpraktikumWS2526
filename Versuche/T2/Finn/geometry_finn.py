# ------------------------------------------------------------
# Geometriedaten und Materialdaten:
# Hier nur direkt gemessene Größen mit Fehlern eintragen
# ------------------------------------------------------------

r0_cfg, dr0_cfg = {}, {}
r_cfg, dr_cfg = {}, {}
F_d_cfg, dF_d_cfg = {}, {}
N_cfg, dN_cfg = {}, {}

# ------------------------------------------------------------
# Abstände Quelle -> Streukörper und Streukörper -> Detektor
# ------------------------------------------------------------

# Ringgeometrie
ring_r0 = {10: 1.30, 20: 0.79, 30: 0.60, 40: 0.39, 50: 0.29}
ring_dr0 = {10: 0.005, 20: 0.005, 30: 0.005, 40: 0.005, 50: 0.005} # Fehler abgeschätzt auf 0.5cm wegen schlechter justage

ring_r = {10: 1.40, 20: 0.57, 30: 0.34, 40: 0.27, 50: 0.22}
ring_dr = {10: 0.005, 20: 0.005, 30: 0.005, 40: 0.005, 50: 0.005}

# Konventionelle Geometrie
conv_r0_val = 0.04955
conv_dr0_val = 0.00003

conv_r_val = 0.0738 + 0.0414 + 0.073
conv_dr_val = np.sqrt(3) * 0.00003

# ------------------------------------------------------------
# Detektorflächen
# ------------------------------------------------------------

# Ringgeometrie:
# F_d = pi * (d_D/2)^2
d_D_ring = 0.0808         # m
dd_D_ring = 0.00003        # m

F_d_ring = np.pi * (d_D_ring / 2)**2
dF_d_ring = F_d_ring * 2 * dd_D_ring / d_D_ring

# Konventionelle Geometrie:
# h = d_L * (a+b+c) / (2*(a+b))
# F_d = pi * h^2

d_L_conv = 0.0245        # m
dd_L_conv = 0.00003      # m

a_conv = 0.0738         # m
da_conv = 0.00003         # m

b_conv = 0.0414          # m
db_conv = 0.00003           # m

c_conv = 0.073           # m
dc_conv = 0.00003           # m

h_conv = d_L_conv * (a_conv + b_conv + c_conv) / (2 * (a_conv + b_conv))

# numerisch oder analytisch; hier analytisch:
dh_ddL = (a_conv + b_conv + c_conv) / (2 * (a_conv + b_conv))
dh_da  = -d_L_conv * c_conv / (2 * (a_conv + b_conv)**2)
dh_db  = -d_L_conv * c_conv / (2 * (a_conv + b_conv)**2)
dh_dc  = d_L_conv / (2 * (a_conv + b_conv))

dh_conv = np.sqrt(
    (dh_ddL * dd_L_conv)**2 +
    (dh_da  * da_conv)**2 +
    (dh_db  * db_conv)**2 +
    (dh_dc  * dc_conv)**2
)

F_d_conv = np.pi * h_conv**2
dF_d_conv = 2 * np.pi * h_conv * dh_conv

# ------------------------------------------------------------
# Elektronenzahlen N
# ------------------------------------------------------------

# Elektronendichten [1/m^3]
n_e_Al = 7.8e29#1.81e29 
dn_e_Al = 0.05e29

n_e_Fe = 2.2e30#1.7e29
dn_e_Fe = 0.2e29 #Quelle????? Unsicherheiten richtig??

# Ring aus Aluminium:
# Beispiel: Ringvolumen V = pi * h * (R_a^2 - R_i^2)

R_a_ring = 0.25         # m
dR_a_ring = 0.01 / np.sqrt(12)        # m

R_i_ring = R_a_ring - 0.0152         # m
dR_i_ring = np.sqrt(0.00003**2 + dR_a_ring**2)        # m

H_ring = 0.015           # m
dH_ring = 0.00003          # m

V_ring = np.pi * H_ring * (R_a_ring**2 - R_i_ring**2)
dV_ring = np.sqrt(
    (np.pi * (R_a_ring**2 - R_i_ring**2) * dH_ring)**2 +
    (np.pi * H_ring * 2 * R_a_ring * dR_a_ring)**2 +
    (np.pi * H_ring * 2 * R_i_ring * dR_i_ring)**2
)

N_ring_Al = n_e_Al * V_ring
dN_ring_Al = np.sqrt(
    (V_ring * dn_e_Al)**2 +
    (n_e_Al * dV_ring)**2
)

# Zylinder für konventionelle Geometrie:
# V = pi * (d/2)^2 * h

# Aluminium-Zylinder
d_cyl_Al = 0.0120
dd_cyl_Al = 0.00003

h_cyl_Al = 0.02105 # Messung 1 wurde genommen
dh_cyl_Al = 0.00003

V_cyl_Al = np.pi * (d_cyl_Al / 2)**2 * h_cyl_Al
dV_cyl_Al = np.sqrt(
    (np.pi * h_cyl_Al * d_cyl_Al / 2 * dd_cyl_Al)**2 +
    (np.pi * (d_cyl_Al / 2)**2 * dh_cyl_Al)**2
)

N_cyl_Al = n_e_Al * V_cyl_Al
dN_cyl_Al = np.sqrt(
    (V_cyl_Al * dn_e_Al)**2 +
    (n_e_Al * dV_cyl_Al)**2
)

# Eisen-Zylinder
d_cyl_Fe = 0.0185
dd_cyl_Fe = 0.00003

h_cyl_Fe = 0.02105
dh_cyl_Fe = 0.00003

V_cyl_Fe = np.pi * (d_cyl_Fe / 2)**2 * h_cyl_Fe
dV_cyl_Fe = np.sqrt(
    (np.pi * h_cyl_Fe * d_cyl_Fe / 2 * dd_cyl_Fe)**2 +
    (np.pi * (d_cyl_Fe / 2)**2 * dh_cyl_Fe)**2
)

N_cyl_Fe = n_e_Fe * V_cyl_Fe
dN_cyl_Fe = np.sqrt(
    (V_cyl_Fe * dn_e_Fe)**2 +
    (n_e_Fe * dV_cyl_Fe)**2
)