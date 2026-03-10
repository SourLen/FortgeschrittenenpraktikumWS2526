import serial
import time
import csv

PORT = "PORT"
BAUD = 115200
OUTFILE = "aufzug.csv"

ser = serial.Serial(PORT, BAUD, timeout=1)
time.sleep(2)
ser.reset_input_buffer()

with open(OUTFILE, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "pi_time",
        "arduino_ms",
        "ax", "ay", "az",
        "mx", "my", "mz"
    ])

    try:
        while True:
            line = ser.readline().decode(errors="ignore").strip()
            if not line:
                continue

            parts = line.split(",")
            if len(parts) != 7:
                continue

            try:
                arduino_ms = int(parts[0])
                ax = float(parts[1])
                ay = float(parts[2])
                az = float(parts[3])
                mx = float(parts[4])
                my = float(parts[5])
                mz = float(parts[6])
            except ValueError:
                continue

            pi_time = time.time()
            writer.writerow([pi_time, arduino_ms, ax, ay, az, mx, my, mz])
            f.flush()

            print(f"{arduino_ms:8d} ms | "
                  f"a=({ax:7.3f}, {ay:7.3f}, {az:7.3f}) m/s² | "
                  f"B=({mx:7.3f}, {my:7.3f}, {mz:7.3f})")
    except KeyboardInterrupt:
        pass

ser.close()