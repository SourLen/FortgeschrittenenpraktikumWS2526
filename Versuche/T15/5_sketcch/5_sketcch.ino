
#include <Adafruit_ADS1X15.h>
#include <Adafruit_MCP4725.h>

Adafruit_MCP4725 dac;
Adafruit_ADS1115 adc;

void setup() 
{
  Serial.begin(9600);
  dac.begin(0x62);

  adc.setGain(GAIN_ONE);
  adc.begin(0x48);

}

void loop() 
{
  int16_t ADC_signal0 = 0;
  int16_t ADC_signal1 = 0;
  int16_t ADC_signal2 = 0;
  int16_t ADC_signal3 = 0;

  uint32_t DAC_signal = 0;


  for(int i=1148; i<= 2784; i = i+4) // gehe die Spannungen von 1.4 bis 3.4 V (v = i * 5 / 4096) in 0.005V = 4bits ( 0.005 * 4096/5) durch
  {
    DAC_signal = i; //12 bit
    dac.setVoltage(DAC_signal, false);
    delay(5);
    ADC_signal0 = adc.readADC_SingleEnded(0);
    ADC_signal1 = adc.readADC_SingleEnded(1);
    ADC_signal2 = adc.readADC_SingleEnded(2);
    ADC_signal3 = adc.readADC_SingleEnded(3);
    Serial.print(i); Serial.print(","); // kommt mit 12 bit an
    Serial.print(ADC_signal0); Serial.print(",");
    Serial.print(ADC_signal1); Serial.print(",");
    Serial.print(ADC_signal2); Serial.print(",");
    Serial.print(ADC_signal3); Serial.print(","); // komm mit 16 bit an
    Serial.println();
    delay(100);
  }
  while (true){}
}