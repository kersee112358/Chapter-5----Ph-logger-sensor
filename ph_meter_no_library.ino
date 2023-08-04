#include <Adafruit_ADS1X15.h>
#include "Nokia_5110.h"
#include <OneWire.h>
#include <DallasTemperature.h>
#include <SPI.h>
#include <SD.h>
#include "RTClib.h"
# define button_pin 5
//Defining the pins for the lcd
#define RST 12
#define CE 11
#define DC 10
#define DIN 9
#define CLK 13
#define VBATPIN A9

Nokia_5110 lcd = Nokia_5110(RST, CE, DC, DIN, CLK);

//32u4 pin for the microSD module
//const int chipSelect = 10;

// Data wire is plugged into digital pin 6 on the Arduino
#define ONE_WIRE_BUS 6

// Setup a oneWire instance to communicate with any OneWire device
OneWire oneWire(ONE_WIRE_BUS);	

// Pass oneWire reference to DallasTemperature library
DallasTemperature sensors(&oneWire);

//Use this for the 16-bit version
Adafruit_ADS1115 ads;  

float voltage,phValue,temperature ;

RTC_PCF8523 rtc;

File myFile;

void setup() {

  sensors.begin();	// Start up the library

  Serial.begin(9600);

  ads.setGain(GAIN_ONE);        // 1x gain   +/- 4.096V  1 bit = 2mV      0.125mV
   if (!ads.begin()) {
    Serial.println("Failed to initialize ADS.");
  }
  
  rtc.begin(); //RTC

  lcd.print("Please Wait ...");
  delay(1000);
  lcd.clear();

  if (!SD.begin(chipSelect)) {
    Serial.println("Card failed, or not present");
    while (1);
  }
  Serial.println("card initialized.");
}

void loop() {

  int16_t adc0;
  float volts0;

  sensors.requestTemperatures(); 
  
  float temperature = sensors.getTempCByIndex(0);

  adc0 = ads.readADC_SingleEnded(0); //reading the voutput of the sensor

  volts0 = ads.computeVolts(adc0); //converting in volts

  lcd.setCursor(0, 0);
  lcd.print("PH: ");
  lcd.println(measure(volts0,temperature));
  lcd.setCursor(0, 1);
  lcd.print("Voltage: ");
  lcd.println(volts0,4);
  lcd.setCursor(0, 2);
  lcd.print("Temp: ");
  lcd.println(sensors.getTempCByIndex(0));
  lcd.setCursor(0, 3);
  lcd.print("Bat: ");
  lcd.print(battery());

  if (digitalRead(button_pin)== HIGH)
    store(measure(volts0,temperature), volts0, temperature);  
delay(1000);
}

float battery(){
  float measuredvbat = analogRead(VBATPIN);
  measuredvbat *= 2;    // we divided by 2, so multiply back
  measuredvbat *= 3.3;  // Multiply by 3.3V, our reference voltage
  measuredvbat /= 1024; // convert to voltage
  Serial.print("VBat: " ); Serial.println(measuredvbat);
  return(measuredvbat);
}

float measure(float volts, float temp){
  float pH=0;

  if (temp <=10)
    pH = -6.27615*volts+16.4456;
  if (temp >10 & temp <=15)
    pH = -6.134969*volts+16.25767;
  if (temp >15 & temp <=20)
    pH = -6.0241*volts+16.09036;
  else
    pH = -6.04351*volts+16.1196615;
  
  Serial.print("Volts: ");
  Serial.println(volts);
  Serial.print("pH :");
  Serial.println(pH);
  return pH;
}

void store(float pH, float volts,float temp){

  myFile = SD.open("pH.txt", FILE_WRITE);
  // if the file opened okay, write to it:
  if (myFile) {

    DateTime now = rtc.now();
    myFile.print(now.year(), DEC);
    myFile.print("/");
    myFile.print(now.month(), DEC);
    myFile.print("/");
    myFile.print(now.day(), DEC);
    myFile.print(" ");
    myFile.print(now.hour(), DEC);
    myFile.print(':');
    myFile.print(now.minute(), DEC);
    myFile.print(':');
    myFile.print(now.second(), DEC); 
    myFile.print(" , ");
    myFile.print(pH,3);
    myFile.print(",");
    myFile.print(volts,3);
    myFile.print(",");
    myFile.print(temp,3);
    myFile.println();
    Serial.println("DONE writing on sd");
    // close the file:
    myFile.close();
  } 
  else {
    // if the file didn't open, print an error:
    Serial.println("error opening test.txt");
  }
}
