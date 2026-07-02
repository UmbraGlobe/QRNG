#include <Arduino.h>
#include <FspTimer.h> 
#include <vector>
#include <analogWave.h>

FspTimer fsp_timer;

const int PIN = A0;                 
const int SINE_FREQ = 287;          

analogWave wave(DAC);


void setup() {
  pinMode(PIN, OUTPUT);
  wave.sine(SINE_FREQ);
}

void loop(){}