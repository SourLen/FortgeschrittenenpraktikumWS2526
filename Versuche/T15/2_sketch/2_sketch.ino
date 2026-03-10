#include <DallasTemperature.h>
#include <OneWire.h>
#define ONE_WIRE_BUS 5
OneWire oneWire(ONE_WIRE_BUS);
DallasTemperature sensors(&oneWire);
int meas_num = 0;
int temp;
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
  serial.print(",");
  serial.println(temp);
  n++;
  delay(1000);
  if (temp < -50 || temp > 100){
    Serial.print("errorrrr")
  }
}
