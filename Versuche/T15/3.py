import serial
import time

ser = serial.Serial('PORT_NAME', 9600)
time.sleep(2)

distances = []
with open('3_output.txt', 'w') as f:
    f.write("Index, deltaT")
    line = ser.readline().decode('utf-8').rstrip()
    if line:
        parts = line.split(',')
        idx = int(parts[0])
        dT = float(parts[1])
        f.write(f"{idx}, {dT}\n") 
ser.close()
        