import serial
import time

ser = serial.Serial('PORT', 9600)
time.sleep(2)

with open('output.txt', 'w') as f:
    count = 0
    max_events = 10
    while count < max_events:
        if ser.in_waiting > 0:
            data = ser.read()
            state = int.from_bytes(data, "little")
            timestampe = time.time()
            f.write(f"{timestampe}: {state}\n")
            count += 1
    ser.close()
    
### Alternative Steuerung delay via python
'''
import serial
import time 

ser = serial.Serial('PORT', 9600)
time.sleep(2)

while True:
    ser.write(b'1')
    time.sleep(1)
    ser.write(b'0')
    time.sleep(1)
'''

