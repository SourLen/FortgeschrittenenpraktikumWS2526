import csv
import matplotlib.pyplot as plt

# Read the data from the CSV file

timestamps = []
mag_x = []
mag_y = []
mag_z = []

accel_x = []
accel_y = []
accel_z = []

with open('noise_data.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        timestamps.append(float(row['Time']))
        mag_x.append(float(row['mag_x']))
        mag_y.append(float(row['mag_y']))
        mag_z.append(float(row['mag_z']))
        accel_x.append(float(row['accel_x']))
        accel_y.append(float(row['accel_y']))
        accel_z.append(float(row['accel_z']))

# 3D Plotting the magnetic data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(mag_x, mag_y, mag_z, label='Magnetic Path', linewidth=0, marker='.')
ax.set_xlabel('Mag X (mT)')
ax.set_ylabel('Mag Y (mT)')
ax.set_zlabel('Mag Z (mT)')
ax.set_title('3D Magnetic Data from LSM9DS1')
ax.legend()
plt.show()

# 3D Plotting the acceleration data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(accel_x, accel_y, accel_z, label='Acceleration Path', linewidth=0, marker='.')
ax.set_xlabel('Accel X (m/s^2)')
ax.set_ylabel('Accel Y (m/s^2)')
ax.set_zlabel('Accel Z (m/s^2)')
ax.set_title('3D Acceleration Data from LSM9DS1')
ax.legend()
plt.show()
