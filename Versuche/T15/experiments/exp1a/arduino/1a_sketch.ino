/*
Authors: Lennart Sauer, Sinan Blanck

Date: Mar. 10, 2026

Control a single LED

*/

int pin = 12;
int ledState = 0;
int baudRate = 9600; //Default
void setup() {
  // put your setup code here, to run once:
  Serial.begin(baudRate);
  pinMode(pin, OUTPUT);
}

void loop() {
  // put your main code here, to run repeatedly:
  ledState = LOW; // Bei Arduino entspricht LOW quasi 0, HIGH quasi 1, total toll
  digitalWrite(pin, ledState);
  Serial.print(millis());
  Serial.print(", ");
  Serial.println(ledState);
  delay(2000);

  ledState = HIGH;
  digitalWrite(pin, ledState);
  Serial.print(millis());
  Serial.print(", ");
  Serial.println(ledState);
  delay(4000);
}
