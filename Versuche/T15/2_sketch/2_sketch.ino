/*
Authors: Lennart Sauer, Sinan Blanck

Date: Mar. 10, 2026

Measure the Temerature with a DS18B20 Sensor

*/

#include <DallasTemperature.h>
#include <OneWire.h>
#define ONE_WIRE_BUS 5
OneWire oneWire(ONE_WIRE_BUS);
DallasTemperature sensors(&oneWire);
int meas_num = 0;
float temp;
void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  sensors.begin();
}

void loop() {
  // put your main code here, to run repeatedly:
  sensors.requestTemperatures();
  temp = sensors.getTempCByIndex(0);

  Serial.print(meas_num);
  Serial.print(",");
  Serial.print(millis());
  Serial.print(",");
  Serial.println(temp);
  meas_num++;
  delay(500);
}