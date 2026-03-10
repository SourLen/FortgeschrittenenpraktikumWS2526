/*
Authors: Lennart Sauer, Sinan Blanck

Date: Mar. 10, 2026

Control two Leds at the same time, letting them blink individually and together

*/

int ledPin1 = 10;
int ledPin2 = 8;

int ledState1;
int ledState2;

void setup() {

  Serial.begin(9600);

  pinMode(ledPin1, OUTPUT);
  pinMode(ledPin2, OUTPUT);

}

void loop() {

  // both OFF
  ledState1 = LOW;
  ledState2 = LOW;

  digitalWrite(ledPin1, ledState1);
  digitalWrite(ledPin2, ledState2);

  Serial.print(millis());
  Serial.print(",");
  Serial.print(ledState1);
  Serial.print(",");
  Serial.println(ledState2);

  delay(1000);

  // LED1 ON
  ledState1 = HIGH;
  ledState2 = LOW;

  digitalWrite(ledPin1, ledState1);
  digitalWrite(ledPin2, ledState2);

  Serial.print(millis());
  Serial.print(",");
  Serial.print(ledState1);
  Serial.print(",");
  Serial.println(ledState2);

  delay(1000);

  // LED2 ON
  ledState1 = HIGH;
  ledState2 = HIGH;

  digitalWrite(ledPin1, ledState1);
  digitalWrite(ledPin2, ledState2);

  Serial.print(millis());
  Serial.print(",");
  Serial.print(ledState1);
  Serial.print(",");
  Serial.println(ledState2);

  delay(1000);

  // both ON
  ledState1 = LOW;
  ledState2 = HIGH;

  digitalWrite(ledPin1, ledState1);
  digitalWrite(ledPin2, ledState2);

  Serial.print(millis());
  Serial.print(",");
  Serial.print(ledState1);
  Serial.print(",");
  Serial.println(ledState2);

  delay(1000);
}