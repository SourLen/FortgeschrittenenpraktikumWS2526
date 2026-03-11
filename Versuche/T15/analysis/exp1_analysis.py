import matplotlib.pyplot as plt
import numpy as np

patha = "../experiments/exp1a/data/1a_output.txt"
pathb = "../experiments/exp1b/data/1b_output.txt"
pathc = "../experiments/exp1c/data/1c_output.txt"

dataa, unita = np.loadtxt(patha, delimiter=",", dtype=float, skiprows=1), "ms"
datab, unitb = np.loadtxt(pathb, delimiter=",", dtype=float, skiprows=1), "ms"
datac, unitc = np.loadtxt(pathc, delimiter=",", dtype=float, skiprows=1), "s"

colors = ["blue", "red", "green"]

def plot_states(data, unit):
    nplot = np.size(data, 1) - 1
    if nplot > 1:
        fig, ax = plt.subplots(nplot+1, 1, sharex=True, figsize = (9,9))
        fig.supxlabel("t in s")
        fig.suptitle(f'Experiment #1{["a","b","c"][nplot-1]}')
        for i in range(nplot):
            if unit == "ms":
                ax[0].step((data[:,0]-data[0,0])/1000, data[:,i+1], alpha = 1/nplot, color = colors[i])
            elif unit == "s":
                ax[0].step((data[:,0]-data[0,0]), data[:,i+1], alpha = 1/nplot, color = colors[i])
            else:
                raise Exception("No valid unit")
            ax[0].set_xlim((0,20))
            ax[0].set_ylabel("State of the LEDs")
        for i in range(nplot):
            if unit == "ms":
                ax[i+1].step((data[:,0]-data[0,0])/1000, data[:,i+1], color = colors[i])
            elif unit == "s":
                ax[i+1].step((data[:,0]-data[0,0]), data[:,i+1], color = colors[i])
            else:
                raise Exception("No valid unit")
            ax[i+1].set_xlim((0,20))
            ax[i+1].set_ylabel(f"State of LED {i+1}")
    else:
        plt.title(f'Experiment #1{["a","b","c"][nplot-1]}')
        if unit == "ms":
            plt.step((data[:,0]-data[0,0])/1000, data[:,1], color = colors[0])
        elif unit == "s":
            plt.step((data[:,0]-data[0,0]), data[:,1], color = colors[0])
        else:
            raise Exception("No valid unit")
        plt.xlim((0,20))
        plt.xlabel("t in s")
        plt.ylabel("State of the LED")
    plt.savefig(f"../temp_figures/exp1/1{["a","b","c"][nplot-1]}.pdf")

plot_states(dataa, unita)
plot_states(datab, unitb)
plot_states(datac, unitc)