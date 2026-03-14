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

def temp_fit(x, T_R, T_0, k):
      return (T_0 - T_R)*np.exp(-k*x) + T_R


# fitting Data 

popt, cov = curve_fit(temp_fit, time, temp, p0=[24.0, 74.1, 0.001], sigma = temp_std)
temp_fit_data = temp_fit(time, popt[0], popt[1], popt[2])

r = temp - temp_fit_data
chisq = np.sum((r/temp_std) ** 2)
ndof = len(temp) - len(popt)

plt.grid()
plt.scatter(time, temp, label = "measurement", marker = ".", color = "black")
plt.plot(time, temp_fit_data, label = f"fit, $\chi^2/ndof =$ {chisq/ndof:.2f}")
plt.title('Temperature measurement and fit')
plt.xlabel('Time in s')
plt.ylabel("Temperature in Kelvin")
plt.legend()
plt.savefig("../temp_figures/exp2/2_1.pdf")
plt.close()

last_30_temp = temp[-130:]
last_30_time = time[-130:]

plt.grid()
plt.scatter(last_30_time, last_30_temp, label = "measurement", marker = ".", color = "black")
plt.title('Temperature measurement of the last 30 measurements')
plt.xlabel('Time in s')
plt.ylabel("Temperature in Kelvin")
plt.legend()
plt.savefig("../temp_figures/exp2/2_2.pdf")
plt.close()

## adding linear Term

def temp_fit_linear(x, T_R, T_0, k, a):
    return (T_0 - T_R)*np.exp(-k*x) + T_R+a*x


popt_linear, cov_linear = curve_fit(temp_fit_linear, time, temp, p0=[24.0, 74.1, 0.001, 1.0], sigma = temp_std)
temp_fit_data_linear = temp_fit_linear(time, popt_linear[0], popt_linear[1], popt_linear[2], popt_linear[3])

r_linear = temp - temp_fit_data_linear
chisq_linear = np.sum((r_linear/temp_std) ** 2)
print(chisq_linear/(ndof-1))

plt.grid()
plt.scatter(time, temp, label = "measurement", marker = ".", color = "black")
plt.plot(time, temp_fit_data_linear, label = f"fit, $\chi^2/ndof =$ {chisq_linear/ndof:.2f}")
plt.title('Temperature measurement and fit')
plt.xlabel('Time in s')
plt.ylabel("Temperature in Kelvin")
plt.legend()
plt.savefig("../temp_figures/exp2/2_3.pdf")
plt.close()

plt.grid()
plt.scatter(time, r_linear, label = "measurement", marker = ".", color = "black")
plt.title('Temperature measurement and fit')
plt.xlabel('Time in s')
plt.ylabel("Temperature in Kelvin")
plt.legend()
plt.savefig("../temp_figures/exp2/2_4.pdf")
plt.close()

print(popt_linear, cov_linear)