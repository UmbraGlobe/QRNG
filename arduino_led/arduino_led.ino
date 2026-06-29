#include <Arduino.h>
#include <FspTimer.h> 
#include <vector>

FspTimer fsp_timer;

const int PIN = A0;                 
const int SINE_FREQ = 500;          
const int SAMPLES_PER_CYCLE = 64;    
const int TIMER_FREQ = SINE_FREQ * SAMPLES_PER_CYCLE; 

const std::vector<int> message = {0,0,0,0,0,1};
int length;

volatile int messageIndex = 0;
volatile int tickCounter = 0;
volatile int sampleIndex = 0;


uint16_t sineTable[SAMPLES_PER_CYCLE];

static void timerCallback(timer_callback_args_t *p_args) {
  tickCounter++;
  if (tickCounter >= TIMER_FREQ) {
    tickCounter = 0;
    messageIndex = (messageIndex + 1) % length;
  }
  
  int code = message[messageIndex];
  
  if (code) {
    analogWrite(PIN, sineTable[sampleIndex]);
    sampleIndex = (sampleIndex + 1) % SAMPLES_PER_CYCLE;
  } else {
    analogWrite(PIN, 0); 
    sampleIndex = 0; 
  }
}

void setup() {
  pinMode(PIN, OUTPUT);
  length = message.size();
  
  analogWriteResolution(12);

  for (int i = 0; i < SAMPLES_PER_CYCLE; i++) {
    float radians = (2.0 * PI * i) / SAMPLES_PER_CYCLE;
    sineTable[i] = (uint16_t)(4095.0 / 2.0 * (1.0 + sin(radians))); 
  }

  uint8_t type = GPT_TIMER;
  int8_t channel = FspTimer::get_available_timer(type);
  if (channel < 0) return; 
  
  fsp_timer.begin(TIMER_MODE_PERIODIC, type, channel, TIMER_FREQ, 50.0, timerCallback, nullptr);
  fsp_timer.setup_overflow_irq(); 
  fsp_timer.open();
  fsp_timer.start();
}

void loop() { 
  Serial.println("going");
}