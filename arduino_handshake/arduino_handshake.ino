#include <Arduino.h>
#include <FspTimer.h> 

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

static void timerCallback(timer_callback_args_t *p_args) {
  fillBuf->packets[sampleIndex].timestamp = micros();
  fillBuf->packets[sampleIndex].channels[0] = analogRead(A0);
  fillBuf->packets[sampleIndex].channels[1] = analogRead(A1);
  fillBuf->packets[sampleIndex].channels[2] = analogRead(A2);
  fillBuf->packets[sampleIndex].channels[3] = analogRead(A3);
  sampleIndex++;
  if (sampleIndex >= PAYLOAD_SAMPLES)
  {
    sendBuf = fillBuf;
    if (fillBuf == &bufA) fillBuf = &bufB;
    else fillBuf = &bufA;
    sampleIndex = 0;
    dataReady = true;
  }
}

void setup() {
  Serial.begin(230400);

  analogReadResolution(14);
  pinMode(LED_BUILTIN, OUTPUT);
  digitalWrite(LED_BUILTIN, LOW);


  //set header bytes: AA BB
  bufA.header = 0xBBAA;
  bufB.header = 0xBBAA;
  uint8_t type = GPT_TIMER;
  int8_t channel = FspTimer::get_available_timer(type);
  if (channel < 0) return; 
  //setup timer
  fsp_timer.begin(TIMER_MODE_PERIODIC, type, channel, 1500, 50.0, timerCallback, nullptr);
  fsp_timer.setup_overflow_irq(); 
  fsp_timer.open();
  fsp_timer.start();
}

void loop() {
  if (dataReady)
  {
    dataReady = false;
    Serial.write((uint8_t*)sendBuf, sizeof(DataPayload));
  }
}