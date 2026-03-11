/*
Author: Lennart Sauer, Sinan Blanck

Date: Mar. 10, 2026

HC-SR04 measurement
VCC = +5V
Trigger = Pin 3
Echo = Pin 2
*/

#define SPEED_OF_SOUND 0.0343  // cm / microsecond

const int echoPin = 2;
const int triggerPin = 3;

unsigned long meas_num = 0;
unsigned long duration_us = 0;
float distance_cm = 0.0;
unsigned long timestamp_ms = 0;

void setup() {
  Serial.begin(9600);
  pinMode(triggerPin, OUTPUT);
  pinMode(echoPin, INPUT);
  digitalWrite(triggerPin, LOW);
}

void loop() {
  // Trigger pulse
  digitalWrite(triggerPin, LOW);
  delayMicroseconds(4);
  digitalWrite(triggerPin, HIGH);
  delayMicroseconds(10);
  digitalWrite(triggerPin, LOW);

  // Echo duration in microseconds
  duration_us = pulseIn(echoPin, HIGH, 30000);  // timeout 30 ms

  // Ignore failed readings
  if (duration_us == 0) {
    return;
  }

  // Distance in cm
  distance_cm = SPEED_OF_SOUND * duration_us / 2.0;

  // Timestamp in ms since Arduino start
  timestamp_ms = millis();

  // One clean CSV line:
  Serial.print(meas_num);
  Serial.print(",");
  Serial.print(duration_us);
  Serial.print(",");
  Serial.print(distance_cm);
  Serial.print(",");
  Serial.println(timestamp_ms);

  meas_num++;

  delay(50);  // important: prevents flooding the serial buffer
}