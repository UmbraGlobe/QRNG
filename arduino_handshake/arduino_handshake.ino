#include <Arduino.h>
#include <FspTimer.h>
#include <analogWave.h>
FspTimer fsp_timer;


volatile bool dataReady = false;

const int PAYLOAD_SAMPLES = 50;

struct __attribute__((packed)) DataPacket {
  uint32_t timestamp;
  uint16_t channels[4];
};


struct __attribute__((packed)) DataPayload {
  uint16_t header;
  DataPacket packets[50];
};

volatile DataPayload bufA;
volatile DataPayload bufB;

volatile DataPayload* volatile fillBuf = &bufA;
volatile DataPayload* volatile sendBuf = nullptr;

unsigned int sampleIndex = 0;

static void timerCallback(timer_callback_args_t* p_args) {
  fillBuf->packets[sampleIndex].timestamp = micros();
  fillBuf->packets[sampleIndex].channels[0] = analogRead(A4);
  fillBuf->packets[sampleIndex].channels[1] = analogRead(A1);
  fillBuf->packets[sampleIndex].channels[2] = analogRead(A2);
  fillBuf->packets[sampleIndex].channels[3] = analogRead(A3);
  sampleIndex++;
  if (sampleIndex >= PAYLOAD_SAMPLES) {
    sampleIndex = 0;
    if (!dataReady) {
      sendBuf = fillBuf;
      if (fillBuf == &bufA) fillBuf = &bufB;
      else fillBuf = &bufA;
      dataReady = true;
    }
  }
}

#define SAMPLE_COUNT 128

const float MAX_VOLTAGE = 1.2;
const float MIN_VOLTAGE = 0.4;
const float SYSTEM_VOLTAGE = 5;

const int DAC_RESOLUTION = 65535;
const int VAL_HIGH = (MAX_VOLTAGE * DAC_RESOLUTION) / SYSTEM_VOLTAGE;
const int VAL_LOW  = (MIN_VOLTAGE * DAC_RESOLUTION) / SYSTEM_VOLTAGE;

const int centerOffset = (VAL_HIGH + VAL_LOW) / 2;
const int amplitude = (VAL_HIGH - VAL_LOW) / 2;
const int SINE_FREQ = 930; 

uint16_t offsetSin[SAMPLE_COUNT];

analogWave wave(DAC, offsetSin, SAMPLE_COUNT, 0);

void setup() {
  Serial.begin(921600);
  pinMode(DAC, OUTPUT);

  for (int i = 0; i < SAMPLE_COUNT; i++) {
    offsetSin[i] = (uint16_t)(centerOffset + (amplitude * sin(2.0 * PI * i / SAMPLE_COUNT)));
  }

  wave.begin(SINE_FREQ);
  analogReadResolution(14);

  //set header bytes: AA BB
  bufA.header = 0xBBAA;
  bufB.header = 0xBBAA;
  uint8_t type = GPT_TIMER;
  int8_t channel = FspTimer::get_available_timer(type);
  if (channel < 0) return;
  //setup timer
  fsp_timer.begin(TIMER_MODE_PERIODIC, type, channel, 2500, 50.0, timerCallback, nullptr);
  fsp_timer.setup_overflow_irq();
  fsp_timer.open();
  fsp_timer.start();
}

int count = 0;

void loop() {
  if (dataReady) {
    count++;
    Serial.write((uint8_t*)sendBuf, sizeof(DataPayload));
    dataReady = false;
    if (count % 100 == 0)
    {
      //wave.freq(SINE_FREQ+count/10);
    }
  }
}