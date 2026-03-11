import serial
import time
import matplotlib.pyplot as plt
import numpy as np

PORT = 'COM5'
BAUD = 9600
OUTPUT_FILE = '3_output_schall_50cm_outside.txt'

ser = serial.Serial(PORT, BAUD, timeout=1)
time.sleep(2)
ser.reset_input_buffer()

with open(OUTPUT_FILE, 'w') as f:
    f.write("Index,deltaT_us,distance_cm,timestamp_ms\n")

    try:
        while True:
            line = ser.readline().decode('utf-8', errors='ignore').strip()

            if not line:
                continue

            parts = line.split(',')
            if len(parts) != 4:
                print("Verworfen:", line)
                continue

            try:
                idx = int(parts[0])
                dT_us = float(parts[1])
                dist_cm = float(parts[2])
                timestamp_ms = float(parts[3])
            except ValueError:
                print("Fehler beim Parsen:", line)
                continue

            f.write(f"{idx},{dT_us},{dist_cm},{timestamp_ms}\n")
            f.flush()

            print(idx, dT_us, dist_cm, timestamp_ms)

    except KeyboardInterrupt:
        print("Messung beendet.")

ser.close()


'''

data = np.loadtxt(OUTPUT_FILE, delimiter=',', skiprows=1)

t = data[:, 3] * 1e-3      # ms -> s
x = data[:, 2] * 1e-2      # cm -> m

# optional: doppelte Zeitwerte entfernen
mask = np.diff(t, prepend=t[0] - 1) > 0
t = t[mask]
x = x[mask]

v = np.gradient(x, t)
a = np.gradient(v, t)

plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(t, x)
plt.title('Distance vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Distance (m)')

plt.subplot(3, 1, 2)
plt.plot(t, v)
plt.title('Velocity vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')

plt.subplot(3, 1, 3)
plt.plot(t, a)
plt.title('Acceleration vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (m/s²)')

plt.tight_layout()
plt.show()


'''
# Daten laden und auswerten
