/*
Author: Lennart Sauer, Sinan Blanck

Date: Mar. 10, 2026

This code writes to the state of an LED, 
which is connected to the Arduino

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
  ledState = LOW; // Bei Arduino entspricht LOW quasi 0, HIGH quasi 1, total toll und logisch haha hihi
  digitalWrite(pin, ledState);
  Serial.write(ledState);
  delay(2000);

  ledState = HIGH;
  digitalWrite(pin, ledState);
  Serial.write(ledState);
  delay(4000);
  /*
  ALternative falls Python die Delays steuern soll
  void loop(){
    if (Serial.available() > 0){
      ledState = Serial.read();

      if (ledState == '1'){
        digitalWrite(pin, HIGH);
      }
      if (ledState == '0'){
        digitalWrite(pin, LOW)
      }

    }

  }
  
  */
}
