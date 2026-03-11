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

with open('circle_data_new.csv', 'r') as f:
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


### Open eleveator_data_new.csv and aufzug_nahe_mb.csv and compare the accel_z in a figure
### make sure to rescale the timestamps appropriately so that they can be compared despite the measurements being taken at different times
### ALSO ADD THE ACCEL_Z DATA FROM THE NEW lastenaufzug.csv FILE TO THE PLOT, RESCALING THE TIMESTAMPS AS WELL

elevator_timestamps = []
elevator_accel_z = []
with open('elevator_data_new.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        elevator_timestamps.append(float(row['Time']))
        elevator_accel_z.append(float(row['accel_z']))
aufzug_timestamps = []
aufzug_accel_z = []
with open('aufzug_nahe_mb.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        aufzug_timestamps.append(float(row['Time']))
        aufzug_accel_z.append(float(row['accel_z']))
        
lastenaufzug_timestamps = []
lastenaufzug_accel_z = []
with open('lastenaufzug.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        lastenaufzug_timestamps.append(float(row['Time']))
        lastenaufzug_accel_z.append(float(row['accel_z']))
        
# Rescale timestamps to start from zero

elevator_timestamps = [t - elevator_timestamps[0] for t in elevator_timestamps]
aufzug_timestamps = [t - aufzug_timestamps[0] for t in aufzug_timestamps]
lastenaufzug_timestamps = [t - lastenaufzug_timestamps[0] for t in lastenaufzug_timestamps]
# Plotting the accel_z data for both datasets
plt.figure()
plt.plot(elevator_timestamps, elevator_accel_z, label='Elevator Data', linewidth=1)
plt.plot(aufzug_timestamps, aufzug_accel_z, label='Aufzug Nahe MB', linewidth=1)
plt.plot(lastenaufzug_timestamps, lastenaufzug_accel_z, label='Lastenaufzug', linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Accel Z (m/s^2)')
plt.title('Comparison of Accel Z between Elevator and Aufzug Nahe MB')
plt.legend()
plt.show()

