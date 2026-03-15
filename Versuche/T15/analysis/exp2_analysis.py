import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat, umath
from scipy.optimize import curve_fit

# Loading the Data
data = np.loadtxt("../experiments/exp2/data/2_temp_output.txt", skiprows = 30, delimiter=",", dtype=float, usecols=[1,2])
noise = np.loadtxt("../experiments/exp2/data/2_noise_output.txt", skiprows = 1, delimiter=",", dtype=float, usecols=2)

time = data[:,0]/1000
temp = data[:,1]
temp_std_manufactorer = 0.5/np.sqrt(3)
temp_std_noise = np.std(noise)
temp_std_dig = 0.065/np.sqrt(12)

temp_std = np.sqrt(temp_std_dig**2 + temp_std_manufactorer**2 + temp_std_dig**2)

# Plotting the Digiatalization

last_100_temp = temp[-100:]
last_100_time = time[-100:]

plt.grid()
plt.scatter(last_100_time, last_100_temp, label = "measurement", marker = ".", color = "black")
plt.title('Temperature measurement of the last 100 measurements')
plt.xlabel('Time in s')
plt.ylabel("Temperature in Kelvin")
plt.legend()
plt.savefig("../temp_figures/exp2/2_1.pdf")
plt.close()





# Fitting Data with the normal model

def temp_fit(x, T_R, T_0, k):
      return (T_0 - T_R)*np.exp(-k*x) + T_R

popt, cov = curve_fit(temp_fit, time, temp, p0=[24.0, 74.1, 0.001], sigma = temp_std)
temp_fit_data = temp_fit(time, popt[0], popt[1], popt[2])

r = temp - temp_fit_data
chisq = np.sum((r/temp_std) ** 2)
ndof = len(temp) - len(popt)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[3,1])
fig.supxlabel("Time in s")
ax1.scatter(time, temp, label = "measurement", marker = ".", color = "black")
ax1.plot(time, temp_fit_data, label = f"fit, $\\chi^2/ndof =$ {chisq/ndof:.2f}")
ax1.set_title('Temperature measurement and fit')
ax1.set_ylabel("Temperature in Kelvin")
ax1.legend()
ax1.grid()
ax2.errorbar(time[::100], r[::100], yerr=temp_std, fmt = ".", color = "black", label="Residues")
ax2.set_title('Residual Plot')
ax2.set_ylabel("Temperature in Kelvin")
ax2.legend()
ax2.grid()
plt.savefig("../temp_figures/exp2/2_2.pdf")
plt.close()



## adding linear Term

def temp_fit_linear(x, T_R, T_0, k, a):
    return (T_0 - T_R*(1+a*x))*np.exp(-k*x) + T_R*(1+a*x)


popt_linear, cov_linear = curve_fit(temp_fit_linear, time, temp, p0=[24.0, 74.1, 0.001, 1.0], sigma = temp_std)
temp_fit_data_linear = temp_fit_linear(time, popt_linear[0], popt_linear[1], popt_linear[2], popt_linear[3])

r_linear = temp - temp_fit_data_linear
chisq_linear = np.sum((r_linear/temp_std) ** 2)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[3,1])
fig.supxlabel("Time in s")
ax1.scatter(time, temp, label = "measurement", marker = ".", color = "black")
ax1.plot(time, temp_fit_data_linear, label = f"fit, $\\chi^2/ndof =$ {chisq_linear/ndof:.2f}")
ax1.set_title('Temperature measurement and fit')
ax1.set_ylabel("Temperature in Kelvin")
ax1.legend()
ax1.grid()
ax2.errorbar(time[::100], r_linear[::100], yerr=temp_std, fmt = ".", color = "black", label="Residues (every 100th)")
ax2.set_title('Residual Plot')
ax2.set_ylabel("Temperature in Kelvin")
ax2.legend()
ax2.grid()
plt.savefig("../temp_figures/exp2/2_3.pdf")
plt.close()


print(popt[0], cov[0,0])
print(popt[1], cov[1,1])
print(popt[2], cov[2,2])