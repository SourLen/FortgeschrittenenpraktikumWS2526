import serial
import time

ser = serial.Serial('COM6', 9600)
time.sleep(2)

with open('1b_output.txt', 'w') as f:
    f.write("Timestamp, State1, State2\n")
    counts = 0
    while counts < 60:
        line = ser.readline().decode('utf-8').rstrip()
        parts = line.split(',')
        timestamp = float(parts[0])
        State1= int(parts[1])
        State2 = int(parts[2])
        f.write(f"{timestamp}, {State1}, {State2}\n")
        counts += 1
ser.close()