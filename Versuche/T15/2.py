import serial
import time

ser = serial.Serial('COM4', 9600)
time.sleep(2)

with open('2_noise_output.txt', 'w') as f:
    f.write("Index, Timestamp, Temperatur in Grad Celsius\n")
    timestamp = 0
    while True:
        line = ser.readline().decode('utf-8').rstrip() 
        parts = line.split(',')
        idx = int(parts[0])
        timestamp = float(parts[1])
        temp = float(parts[2])
        print(f"Index: {idx}, Timestamp: {timestamp}, Temperature: {temp}")
        f.write(f"{idx}, {timestamp}, {temp}\n")
ser.close()