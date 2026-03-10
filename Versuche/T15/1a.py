import serial
import time

ser = serial.Serial('COM4', 9600)
time.sleep(2)

with open('1a_output.txt', 'w') as f:
    f.write("Timestamp, State\n")
    counts = 0
    while counts < 10:
        line = ser.readline().decode('utf-8').rstrip()
        parts = line.split(',')
        timestamp = float(parts[0])
        State = int(parts[1])
        f.write(f"{timestamp}, {State}\n")
        counts += 1
ser.close()

