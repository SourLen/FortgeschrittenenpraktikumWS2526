import time
import board
import busio
import adafruit_lsm9ds1
import RPi.GPIO as GPIO
import csv

# Define the switch
switch = 14   # GPIO14 = physical pin 8

GPIO.setmode(GPIO.BCM)
GPIO.setup(switch, GPIO.IN, pull_up_down=GPIO.PUD_UP)

# I2C setup with LSM9DS1
i2c = busio.I2C(board.SCL, board.SDA)
sensor = adafruit_lsm9ds1.LSM9DS1_I2C(i2c)

# Setting up csv

f = open("elevator_data.csv", "w")

writer = csv.writer(f)
writer.writerow(["Time", "accel_x", "accel_y", "accel_z", "gyro_x", "gyro_y", "gyro_z", "mag_x", "mag_y", "mag_z"])


print("Waiting for switch pin to go LOW...")

try:
	while True:
		if GPIO.input(switch) == GPIO.LOW:  # pin down
			accel_x, accel_y, accel_z = sensor.acceleration
			gyro_x, gyro_y, gyro_z = sensor.gyro
			mag_x, mag_y, mag_z = sensor.magnetic
			timestamp = time.time() 

			print("Acceleration (m/s^2):", accel_x, accel_y, accel_z)
			print("Gyroscope (rad/s):", gyro_x, gyro_y, gyro_z)
			print("Magnetometer (gauss):", mag_x, mag_y, mag_z)
			print("------------------")

			row = [timestamp, accel_x, accel_y, accel_z, gyro_x, gyro_y, gyro_z,mag_x, mag_y, mag_z]
			writer.writerow(row)
		else:
			print("Pin HIGH -> not measuring")

		time.sleep(0.5)

except KeyboardInterrupt:
	GPIO.cleanup()
	f.close()