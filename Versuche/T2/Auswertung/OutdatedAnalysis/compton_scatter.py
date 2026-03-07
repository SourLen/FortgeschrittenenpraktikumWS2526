from praktikum import analyse
import numpy as np
import matplotlib.pyplot as plt

ring_noise_5_10deg = np.loadtxt("../Messdaten/ComptonScattering/Cs137_41B_ring_10_degrees_5min_noise.tka")[2:]
ring_30_10deg= np.loadtxt("../Messdaten/ComptonScattering/Cs137_41B_ring_10_degrees_30min.tka")[2:]
ring_20_20deg = np.loadtxt("../Messdaten/ComptonScattering/Cs137_41B_ring_20_degrees_20min.tka")[2:]
ring_noise_5_20deg = np.loadtxt("../Messdaten/ComptonScattering/Cs137_41B_ring_20_degrees_5min_noise.tka")[2:]
x_ring_noise_5_10deg = np.arange(len(ring_noise_5_10deg))
x_ring_30_10deg = np.arange(len(ring_30_10deg))
x_ring_20_20deg = np.arange(len(ring_20_20deg))
x_ring_noise_5_20deg = np.arange(len(ring_noise_5_20deg))
ring_30_10deg_no_noise = ring_30_10deg - 6*ring_noise_5_10deg
ring_20_20deg_no_noise = ring_20_20deg - 4*ring_noise_5_20deg
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15), sharex=True)
ax1.plot(x_ring_noise_5_10deg, ring_noise_5_10deg, label="Noise", marker=".")
ax1.plot(x_ring_noise_5_20deg, ring_noise_5_20deg, label="Noise", marker=".")
ax2.plot(x_ring_30_10deg, ring_30_10deg, label="Ring 10 degrees 30 min", marker=".")
ax2.plot(x_ring_20_20deg, ring_20_20deg, label="Ring 20 degrees 20 min", marker=".")
ax2.set_ylabel("counts")
ax3.plot(x_ring_20_20deg, ring_20_20deg_no_noise, label="Ring 20 degrees 20 min no noise", marker=".")
ax1.set_ylabel("counts")
ax2.set_ylabel("counts")
ax3.set_ylabel("counts")
ax3.set_xlabel("channel")
ax1.legend()
ax2.legend()
ax3.legend()

plt.tight_layout()
plt.show()