# QRNG analysis scripts (split from main.ipynb)

Each notebook section is now its own script. The code is copied from the notebook cells.
Run a script from this folder, for example `python 06_single_sensor_stability.py`.

## Files

| File | Notebook cell | What it does |
|---|---|---|
| `qrng_common.py` | Global Setup | Config (PORT, BAUD, LOAD_FROM_FILE, file paths, sensors), serial and file readers, save helpers. **Edit the config here.** |
| `qrng_viz.py` | Visualization | Plot styling and summary figures |
| `01_calibrate.py` | Calibrate | Measures Fs and carrier, saves `calibration.json` |
| `02_capture_raw.py` | Capture and Save | Live capture split into CSV trials |
| `03_live_waveform.py` | Live Waveform | Live window |
| `04_spectrogram.py` | Spectrogram | |
| `05_live_psd.py` | Live Power Spectrum | Live window |
| `06_single_sensor_stability.py` | Single-Sensor Stability | |
| `07_shared_led_noise.py` | Shared LED Fluctuation | |
| `08_rls_common_mode.py` | Online Common-Mode (RLS) | Needs `padasip` |
| `09_fingerprint_averaging.py` | Fingerprint Averaging | |
| `10_live_pairwise_bits.py` | Live Pairwise Bits | Live window |
| `11_aliasing_check.py` | Aliasing Check | |
| `12_noise_floor_snr.py` | Noise Floor and SNR | |
| `13_envelope_ripple.py` | Envelope Ripple | Live window |
| `14_dark_test.py` | Dark Test | |

Run `01_calibrate.py` first (and again when the data file or sample rate changes). The other scripts read
`calibration.json` for Fs and CARRIER_FREQ. Without it they print a warning and use Fs = 1500 Hz,
CARRIER_FREQ = 0.

## Changes from the notebook code

1. Shared values: Fs and CARRIER_FREQ go through `calibration.json` instead of the notebook kernel (added lines at the end of `qrng_common.py` and of `01_calibrate.py`).
2. In `03`, `05`, `08`, `09`, `10`, `13`, the line `with serial_port():` became `with serial_port() as ser:`. The notebook shared `ser` between cells. A script needs it bound locally, and the setup cell's own comment recommends this form. Only live (serial) mode is affected.
3. `qrng_viz.py` ends with `from qrng_common import VIZ_DPI`. In the notebook the setup cell ran after the plotting cell and set VIZ_DPI = 160, so this keeps that value.
4. Each script starts with its notebook description as a docstring plus the imports it needs.

Nothing else in the code was changed.

## Before running on the laptop

- `qrng_common.py` has `PORT = "COM8"` and `Path(r"paper_data\data.csv")`. Both are Windows-style. On macOS or Linux the port is something like `/dev/cu.usbmodem...` or `/dev/ttyACM0`, and the path needs forward slashes.
- The path is relative, so run scripts from the folder that contains `paper_data`.
- `import serial` needs `pip install pyserial`.
- The cells set `CARRIER_FREQ = 930.5` by hand after estimating it (in `01_calibrate.py` and again in `14_dark_test.py`).
