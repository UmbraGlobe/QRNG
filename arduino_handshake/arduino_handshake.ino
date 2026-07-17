#include <Arduino.h>
#include <FspTimer.h>
#include <analogWave.h>
#include "analog.h"  // Added for direct Renesas ADC access
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

static adc_instance_ctrl_t adcCtrl{};
static adc_extended_cfg_t adcExt{};
static adc_cfg_t adcCfg{};
static adc_channel_cfg_t adcChannelCfg{};

static const uint8_t ADC_CH_A1 = 0;
static const uint8_t ADC_CH_A2 = 1;
static const uint8_t ADC_CH_A3 = 2;
static const uint8_t ADC_CH_A4 = 21;

static bool initializeAdc() {

  // Put Arduino pins A1–A4 into analog mode.
  pinPeripheral(
    digitalPinToBspPin(A1),
    IOPORT_CFG_ANALOG_ENABLE
  );

  pinPeripheral(
    digitalPinToBspPin(A2),
    IOPORT_CFG_ANALOG_ENABLE
  );

  pinPeripheral(
    digitalPinToBspPin(A3),
    IOPORT_CFG_ANALOG_ENABLE
  );

  pinPeripheral(
    digitalPinToBspPin(A4),
    IOPORT_CFG_ANALOG_ENABLE
  );

  // Disable hardware averaging.
  adcExt.add_average_count = ADC_ADD_OFF;

  // Do not erase a result register after it is read.
  adcExt.clearing = ADC_CLEAR_AFTER_READ_OFF;

  // Group B is not used, but this field needs a valid setting.
  adcExt.trigger_group_b = ADC_TRIGGER_SYNC_ELC;

  // Disable double-trigger mode.
  adcExt.double_trigger_mode =
    ADC_DOUBLE_TRIGGER_DISABLED;

  // Use AVCC0 and AVSS0 as the ADC voltage references.
  adcExt.adc_vref_control =
    ADC_VREF_CONTROL_AVCC0_AVSS0;

  // Use the normal per-channel ADC result registers.
  adcExt.enable_adbuf = 0;

  // Disable ADC comparison-window interrupts.
  adcExt.window_a_irq = FSP_INVALID_VECTOR;
  adcExt.window_a_ipl = 12;

  adcExt.window_b_irq = FSP_INVALID_VECTOR;
  adcExt.window_b_ipl = 12;

  adcCfg.unit = 0;

  // Each trigger converts all four selected channels once.
  adcCfg.mode = ADC_MODE_SINGLE_SCAN;

  adcCfg.resolution = ADC_RESOLUTION_14_BIT;

  adcCfg.alignment =
    static_cast<adc_alignment_t>(
      ADC_ALIGNMENT_RIGHT
    );

  // readAdcFrame() starts each scan through software.
  adcCfg.trigger = ADC_TRIGGER_SOFTWARE;

  // No ADC scan-complete callback.
  adcCfg.p_callback = nullptr;
  adcCfg.p_context = nullptr;

  adcCfg.p_extend = &adcExt;

  // Scan completion is checked through polling.
  adcCfg.scan_end_irq = FSP_INVALID_VECTOR;
  adcCfg.scan_end_ipl = 12;

  adcCfg.scan_end_b_irq = FSP_INVALID_VECTOR;
  adcCfg.scan_end_b_ipl = 12;

  adcChannelCfg.scan_mask =
    (1UL << ADC_CH_A1) |
    (1UL << ADC_CH_A2) |
    (1UL << ADC_CH_A3) |
    (1UL << ADC_CH_A4);

  adcChannelCfg.scan_mask_group_b = 0;

  // No hardware averaging.
  adcChannelCfg.add_mask = 0;

  // Group scanning is not used.
  adcChannelCfg.priority_group_a =
    ADC_GROUP_A_PRIORITY_OFF;

  // Dedicated sample-and-hold is not used.
  adcChannelCfg.sample_hold_mask = 0;
  adcChannelCfg.sample_hold_states = 0;

  // Window comparison is not used.
  adcChannelCfg.p_window_cfg = nullptr;

  if (R_ADC_Open(&adcCtrl, &adcCfg) != FSP_SUCCESS) {
    return false;
  }

  if (R_ADC_ScanCfg(
        &adcCtrl,
        &adcChannelCfg
      ) != FSP_SUCCESS) {

    R_ADC_Close(&adcCtrl);
    return false;
  }

  return true;
}

static bool readAdcFrame(uint16_t output[4]) {

  // Start one scan containing A1, A2, A3 and A4.
  if (R_ADC_ScanStart(&adcCtrl) != FSP_SUCCESS) {
    return false;
  }

  adc_status_t status{};

  // Prevent an ADC problem from trapping the timer interrupt forever.
  uint32_t timeoutGuard = 10000;

  do {
    if (R_ADC_StatusGet(
          &adcCtrl,
          &status
        ) != FSP_SUCCESS) {
      return false;
    }

    if (timeoutGuard == 0) {
      return false;
    }

    timeoutGuard--;

  } while (
    status.state == ADC_STATE_SCAN_IN_PROGRESS
  );

  if (R_ADC_Read(
        &adcCtrl,
        static_cast<adc_channel_t>(ADC_CH_A4),
        &output[0]
      ) != FSP_SUCCESS) {
    return false;
  }

  if (R_ADC_Read(
        &adcCtrl,
        static_cast<adc_channel_t>(ADC_CH_A1),
        &output[1]
      ) != FSP_SUCCESS) {
    return false;
  }

  if (R_ADC_Read(
        &adcCtrl,
        static_cast<adc_channel_t>(ADC_CH_A2),
        &output[2]
      ) != FSP_SUCCESS) {
    return false;
  }

  if (R_ADC_Read(
        &adcCtrl,
        static_cast<adc_channel_t>(ADC_CH_A3),
        &output[3]
      ) != FSP_SUCCESS) {
    return false;
  }

  return true;
}



static void timerCallback(timer_callback_args_t* p_args) {
  fillBuf->packets[sampleIndex].timestamp = micros();

  uint16_t adcValues[4];

    if (readAdcFrame(adcValues)) {

    fillBuf->packets[sampleIndex].channels[0] =
      adcValues[0];  // A4

    fillBuf->packets[sampleIndex].channels[1] =
      adcValues[1];  // A1

    fillBuf->packets[sampleIndex].channels[2] =
      adcValues[2];  // A2

    fillBuf->packets[sampleIndex].channels[3] =
      adcValues[3];  // A3

  } else {
    fillBuf->packets[sampleIndex].channels[0] = 0xFFFF;
    fillBuf->packets[sampleIndex].channels[1] = 0xFFFF;
    fillBuf->packets[sampleIndex].channels[2] = 0xFFFF;
    fillBuf->packets[sampleIndex].channels[3] = 0xFFFF;
  }
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

const float MAX_VOLTAGE = 1.4;
const float MIN_VOLTAGE = 0.6;
const int SYSTEM_VOLTAGE = 5;

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
  
  if (!initializeAdc()) {
    return;
  }
  //set header bytes: AA BB
  bufA.header = 0xBBAA;
  bufB.header = 0xBBAA;
  uint8_t type = GPT_TIMER;
  int8_t channel = FspTimer::get_available_timer(type);
  if (channel < 0) return;
  //setup timer
  fsp_timer.begin(TIMER_MODE_PERIODIC, type, channel, 5000, 50.0, timerCallback, nullptr);
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