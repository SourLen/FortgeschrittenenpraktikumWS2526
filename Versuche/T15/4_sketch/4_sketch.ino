/*
Author,,,,

Date,...

Beschleunigungsmessung, Aufzug und Magnetfeld


*/


#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_LSM9DS1.h>

Adafruit_LSM9DS1 lsm = Adafruit_LSM9DS1();

void setup() {
  Serial.begin(115200);

  if (!lsm.begin()) {
    while (1) {
      delay(10);
    }
  }
  lsm.setupAccel(lsm.LSM9DS1_ACCELRANGE_2G);
  lsm.setupMag(lsm.LSM9DS1_MAGGAIN_4GAUSS);
  lsm.setupGyro(lsm.LSM9DS1_GYROSCALE_245DPS);
}

void loop() {
  lsm.read();

  sensors_event_t accel, mag, gyro, temp;
  lsm.getEvent(&accel, &mag, &gyro, &temp);

  unsigned long t_ms = millis();

  // CSV: t_ms,ax,ay,az,mx,my,mz
  Serial.print(t_ms);
  Serial.print(",");
  Serial.print(accel.acceleration.x, 6);
  Serial.print(",");
  Serial.print(accel.acceleration.y, 6);
  Serial.print(",");
  Serial.print(accel.acceleration.z, 6);
  Serial.print(",");
  Serial.print(mag.magnetic.x, 6);
  Serial.print(",");
  Serial.print(mag.magnetic.y, 6);
  Serial.print(",");
  Serial.println(mag.magnetic.z, 6);

  delay(50);   