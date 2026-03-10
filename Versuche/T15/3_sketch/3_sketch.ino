/*
Author: Lennart Sauer, Sinan Blanck

Date: Mar. 10, 2026

Measuring the speed of Sound using, HC-SR04 Sensor, CVV = +5VDC, Trigger = Pin3, Echo Pin2

*/

int echo = 2;
int trigger = 3;
int meas_num = 0;
long time;
float distance;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  pinMode(trigger, OUTPUT);
  pinMode(echo, INPUT);

}

void loop() {
  // put your main code here, to run repeatedly:
  digitalWrite(trigger, LOW);
  delayMicroseconds(4);
  digitalWrite(trigger, HIGH);
  delayMicroseconds(10);
  digitalWrite(trigger, LOW);

  time = pulseIn(echo, HIGH);
  Serial.println(time)
  meas_num ++;

}
