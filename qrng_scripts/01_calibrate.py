r"""
## Calibrate Sample Rate and Carrier

**Purpose:** Reads a short record, measures the sample rate from its timestamps, checks the dominant tone on each channel, and sets the shared `Fs` and `CARRIER_FREQ` used by later sections.

**Usage:** With the optical carrier on, choose `SAMPLES_TO_CAPTURE` and run this cell after Global Setup. Live mode captures fresh samples; saved-file mode uses the first path in `LOAD_FILE_PATHS`. Confirm that the printed values agree with the hardware before continuing. The current notebook intentionally pins `CARRIER_FREQ` to **930.5 Hz** after checking the data.

Run this again after changing hardware, firmware timing, carrier, or input files. **Do not calibrate from a dark/LED-off record**: without the optical carrier, the peak detector can lock onto a packet or switching spur. The Dark Test pins its own carrier frequency for this reason.

Source: notebook cell 5 (code copied unchanged unless noted below).
"""

from qrng_common import *

### CONFIG

SAMPLES_TO_CAPTURE = 20000  # total number of samples to capture

###

data = {}
Fs = 1500.0  # initial guess for sampling rate
CARRIER_FREQ = 0  # initial guess for carrier frequency in Hz before estimation

data = read_from_file_instant(LOAD_FILE_PATHS[0]) if LOAD_FROM_FILE else read_data_realtime(SAMPLES_TO_CAPTURE)
data = data[:SAMPLES_TO_CAPTURE]  # limit to the first N samples

if len(data) < 2:
    raise RuntimeError("At least two samples are required to estimate the sampling rate.")
dt_us = (np.diff(data['timestamp'].astype(np.uint64)) + (1 << 32)) % (1 << 32)
if not np.all(dt_us > 0):
    raise ValueError("Timestamps must advance strictly between samples.")
Fs = float(1.0 / (np.median(dt_us) * 1e-6))
estimated_carrier_freqs = {}
for s in SENSORS:
    x = adc_to_mv(data[s])
    fft_vals = fft.fft(x)
    f = fft.fftfreq(len(x), d=1/Fs)
    mags = 2 * np.abs(fft_vals) / len(x)
    inband = (f > 14)  # skip the DC / drift / 1-f region
    peaks, _ = signal.find_peaks(mags[inband], distance=5, height=0.1)  # adjust height threshold as needed
    sorted_peaks = sorted(peaks, key=lambda p: mags[inband][p], reverse=True)
    estimated_carrier_freqs[s] = f[inband][sorted_peaks[0]] if sorted_peaks else 0.0

CARRIER_FREQ = np.mean(list(estimated_carrier_freqs.values()))
CARRIER_FREQ = 930.5
print(f"Estimated sampling rate: {Fs:.2f} Hz")
print(f"Average estimated carrier frequency across sensors: {CARRIER_FREQ:.6f} Hz")

# --- added when splitting the notebook ---
# Save the results so the other scripts can read them.
CALIBRATION_FILE.write_text(json.dumps({"Fs": Fs, "CARRIER_FREQ": float(CARRIER_FREQ)}, indent=2))
print(f"Saved calibration to {CALIBRATION_FILE}")
