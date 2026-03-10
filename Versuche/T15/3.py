import serial
import time

ser = serial.Serial('PORT_NAME', 9600)
time.sleep(2)

distances = []
times = []
timestamps = []
with open('3_output.txt', 'w') as f:
    f.write("Timestamp, Distanz, deltaT, Index")
    line = ser.readline().decode('utf-8').rstrip()
    if line:
        parts = line.split(',')
        idx = int(parts[0])
        dist = float(parts[1])
        dT = float(parts[2])
        timestamp = time.time()
        f.write(f"{timestamp}, {dist}, {dT}, {idx}\n") 
    ser.close()
        