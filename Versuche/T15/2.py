import serial
import time

ser = serial.Serial('PORT_NAME', 9600)
time.sleep(2)

with open('2_temp_output.txt', 'w') as f:
    f.write("Timestamp, Temperatur in Grad Celsius, Index\n")
    while True:
        line = ser.readline().decode('utf-8').rstrip() 
        parts = line.split(',')
        idx = int(parts[0])
        temp = float(parts[1])
        timestamp = time.time()
        f.write(f"{timestamp}, {temp}, {idx}\n")
ser.close()