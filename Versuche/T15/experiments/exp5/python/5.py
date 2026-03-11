import serial
import csv
import time
import numpy as np
import matplotlib.pyplot as plt

COLOR = "RED"
PORT = "COM4"          
BAUD = 9600
OUTPUT_FILE = f"./T15/5_measurements/raw_data_{COLOR}.csv"

ADC_RANGE = 4.096      # Volt bei GAIN_ONE
ADC_RES = 32768        # 16 bit

DAC_REF = 5.0
DAC_RES = 4096

R = 1000  

I_threshold = 1e-4

# =============================
# Serielle Verbindung
# =============================

print("Opening serial connection...")
ser = serial.Serial(PORT, BAUD)
time.sleep(2)

data = []

print("Recording data...")

while True:
    line = ser.readline().decode().strip()

    if not line:
        continue

    try:
        parts = line.split(",")

        dac = int(parts[0])
        adc0 = int(parts[1])
        adc1 = int(parts[2])
        adc2 = int(parts[3])
        adc3 = int(parts[4])

        data.append([dac, adc0, adc1, adc2, adc3])

        print(parts)

        if dac >= 2784:
            break

    except:
        pass

ser.close()

print("Measurement finished")

# =============================
# CSV speichern
# =============================

with open(OUTPUT_FILE, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["DAC", "ADC0", "ADC1", "ADC2", "ADC3"])
    writer.writerows(data)

print("Raw data saved to", OUTPUT_FILE)



data = np.array(data)

dac = data[:,0]
adc0 = data[:,1]
adc1 = data[:,2]
adc2 = data[:,3]
adc3 = data[:,4]

V_dac = dac / DAC_RES * DAC_REF

V0 = adc0 / ADC_RES * ADC_RANGE
V1 = adc1 / ADC_RES * ADC_RANGE
V2 = adc2 / ADC_RES * ADC_RANGE
V3 = adc3 / ADC_RES * ADC_RANGE

# LED Spannung
V_led = V0 - V1

# Strom durch Widerstand
I = (V2 - V3) / R

# =============================
# Schwelle bestimmen
# =============================

idx = np.argmin(np.abs(I - I_threshold))
V_threshold = V_led[idx]

print("Threshold voltage (I=1e-4 A):", V_threshold, "V")


plt.figure()

plt.plot(V_led, I, ".")
plt.axhline(I_threshold)
plt.axvline(V_threshold)

plt.xlabel("LED Voltage (V)")
plt.ylabel("Current (A)")
plt.title("LED I-V curve")
#plt.close()
### Save plot to 5_measurements folder
plt.savefig(f"./T15/5_measurements/iv_curve_{COLOR}.png", dpi=300)
plt.show()