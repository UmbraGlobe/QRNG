r"""
## Global Setup: Choose Live Capture or Saved Data

**Purpose:** Loads the libraries and defines the data format, sensors, serial connection, file loading, unit conversion, and save functions used by every later section.

**Choose one input mode before running this cell:**

- **Live hardware:** set `LOAD_FROM_FILE = False`, connect the Arduino, and set `PORT` to its serial port. Keep `BAUD = 921600` unless the firmware uses another rate.
- **Saved data:** set `LOAD_FROM_FILE = True` and place one or more CSV paths in `LOAD_FILE_PATHS`. A compatible raw file contains `timestamp`, `A2V-16_A1`, `A2V-16_A2`, `A2V-16_A3`, and `A2V-16_A4`.

Leave `SENSORS`, `ADC_MAX`, `V_REF`, and `SAMPLES_PER_PACKET` unchanged unless the hardware or firmware format changed. Run Setup again whenever you switch ports, input modes, sensors, or files, then run **Calibrate Sample Rate and Carrier**. Run only one live/interactive section at a time so cells do not compete for the serial port.

Live batch analyses save results under `data/<analysis>/<date>/run_<sample-count>_samples_<time>/`. Saved-file analyses print and display results but do not create another run folder. See [ANALYSIS_SUITE_GUIDE.md](ANALYSIS_SUITE_GUIDE.md) for the complete first-time workflow.

Source: notebook cell 3 (code copied unchanged unless noted below).
"""

import json

import serial
import time
import contextlib
from datetime import datetime
from pathlib import Path
from scipy.signal import get_window
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import numpy.fft as fft
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore
from scipy.signal import welch
import padasip as pa
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, to_rgb


### CONFIG

#PORT = "COM8"
PORT = "/dev/cu.usbmodem3CDC75F0103C2" # this line if for MAC only 
BAUD = 921600
#Ensure all data files have the same sampling rate and carrier freq.

LOAD_FROM_FILE = False
LOAD_FILE_PATHS= [
    Path(__file__).resolve().parent / "paper_data" / "data.csv",
    ]

SENSORS = ['A2V-16_A1', 'A2V-16_A2', 'A2V-16_A3', 'A2V-16_A4']
ADC_MAX = 16383 #2^14 - 1
V_REF = 5
SAMPLES_PER_PACKET = 50

### SUMMARY VISUALIZATION ###

VIZ_DARK = True         
VIZ_DPI  = 160     

###

if not LOAD_FROM_FILE:
    LOAD_FILE_PATHS = []

payload_dtype = np.dtype([
    ('timestamp', '<u4'),
    ('A2V-16_A1', '<u2'),
    ('A2V-16_A2', '<u2'),
    ('A2V-16_A3', '<u2'),
    ('A2V-16_A4', '<u2')
])

BYTES_PER_SAMPLE = payload_dtype.itemsize * SAMPLES_PER_PACKET
ser = None
reader = None

def fir_lowpass_filter(data, cutoff_hz, fs, num_taps=101):
    b = signal.firwin(num_taps, cutoff=cutoff_hz, fs=fs, window='hamming')
    filtered_data = np.asarray(signal.lfilter(b, 1.0, data))
    return filtered_data, b

def adc_to_mv(adc_value):
    return (adc_value / ADC_MAX) * V_REF * 1000  # Convert to millivolts

def fmt_matrix(M, title, fmt="{:>12.4f}"):
    s = f"\n{title}\n"
    s += " " * 10 + "".join(f"{ch:>12}" for ch in SENSORS) + "\n"
    for i, ch in enumerate(SENSORS):
        s += f"{ch:>10}" + "".join(fmt.format(M[i, j]) for j in range(len(SENSORS))) + "\n"
    return s

# reads one 50-sample packet from serial.
# Fix 3: resync by scanning to the next 0xAA 0xBB frame boundary with read_until, so a
# misaligned/corrupt stream self-realigns instead of blindly jumping a fixed 602 bytes
# (which could step past the real boundary forever on false-locks inside the payload).
def try_read_realtime(ser):
    sync = ser.read_until(b'\xAA\xBB')          # discard junk up to the next frame boundary
    if not sync.endswith(b'\xAA\xBB'):
        return None                             # timed out before a sync word appeared
    raw_data = ser.read(BYTES_PER_SAMPLE)
    if len(raw_data) != BYTES_PER_SAMPLE:
        return None
    data = np.frombuffer(raw_data, dtype=payload_dtype)
    dts = np.diff(data['timestamp'].astype(np.int64))
    ok = (all((data[s] <= ADC_MAX).all() for s in SENSORS) and
          (dts > 0).all() and (dts < 5000).all())    # rejects false locks inside payload data
    return data if ok else None

# Read one 50-sample packet from file. With pace=True, preserve the recorded timing
# for GUI diagnostics. reader.exhausted distinguishes EOF from 'not due yet'.
def read_from_file_realtime(file_path, pace=True):
    data = read_from_file_instant(file_path)
    n_packets = len(data) // SAMPLES_PER_PACKET
    dropped = len(data) % SAMPLES_PER_PACKET
    if dropped:
        print(f"Replay ignores {dropped} trailing sample(s) that do not form a full packet.")

    # Accumulated unsigned timestamp deltas remain valid across a uint32 wrap.
    ts = data['timestamp'].astype(np.uint64)
    if len(ts):
        dt = (np.diff(ts) + (1 << 32)) % (1 << 32)
        due_us = np.concatenate(([0.0], np.cumsum(dt, dtype=float)))[::SAMPLES_PER_PACKET]
    else:
        due_us = np.empty(0, dtype=float)

    idx = 0
    wall_start = None

    def reader(*_):
        nonlocal idx, wall_start
        if idx >= n_packets:
            reader.exhausted = True
            return None
        now = time.monotonic()
        if wall_start is None:
            wall_start = now
        if pace and (now - wall_start) * 1e6 < due_us[idx]:
            return None
        start = idx * SAMPLES_PER_PACKET
        idx += 1
        return data[start:start + SAMPLES_PER_PACKET]

    reader.exhausted = (n_packets == 0)
    reader.total_packets = n_packets
    return reader


def read_from_file_instant(file_path):
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Data file not found: {file_path.resolve()}"
        )

    df = pd.read_csv(file_path)

    required_columns = list(payload_dtype.names or ())
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(
            f"{file_path} is missing required columns: {missing}"
        )

    values = df[required_columns].apply(pd.to_numeric, errors='coerce')
    if values.isna().any().any():
        raise ValueError(f"{file_path} contains missing or non-numeric sample values.")
    ts = values['timestamp'].to_numpy()
    if ((ts < 0) | (ts > np.iinfo(np.uint32).max) | (ts != np.floor(ts))).any():
        raise ValueError(f"{file_path} contains timestamps outside the uint32 range.")
    for s in SENSORS:
        x = values[s].to_numpy()
        if ((x < 0) | (x > ADC_MAX) | (x != np.floor(x))).any():
            raise ValueError(f"{file_path} contains invalid ADC values in {s}.")

    return values.to_records(index=False).astype(payload_dtype)

def get_save_dir(base_folder: str, sample_count: int) -> Path:
    now = datetime.now()
    date_str = now.strftime('%m_%d_%Y') 
    time_str = now.strftime('%H_%M_%S')
    save_dir = Path("data") / base_folder / date_str / f"run_{sample_count}_samples_{time_str}"
    save_dir.mkdir(parents=True, exist_ok=True)
    return save_dir

def try_initialize_serial():
    if not LOAD_FROM_FILE:
        global ser
        if ser is not None and ser.is_open:
            ser.close()
        ser = serial.Serial(PORT, BAUD, timeout=1)
        ser.reset_input_buffer()    # Fix 2: drop the stale free-running backlog on open
        ser.reset_output_buffer()

def try_close_serial():
    global ser
    if ser is not None and ser.is_open:
        ser.close()

def drain_to_latest(ser):
    # Fix 2: discard any buffered backlog so capture starts on fresh, frame-aligned data.
    if ser is not None:
        ser.reset_input_buffer()

@contextlib.contextmanager
def serial_port():
    # Fix 1: guaranteed-close port. Use `with serial_port() as ser:` so the handle is
    # released even on KeyboardInterrupt, a thrown exception, or a closed plot window --
    # no more stuck COM port that needs an unplug/replug. Sets the global `ser` so any
    # code that still references it keeps working. Yields None when replaying from file.
    global ser
    if LOAD_FROM_FILE:
        yield None
        return
    ser = serial.Serial(PORT, BAUD, timeout=1)
    try:
        ser.reset_input_buffer()    # Fix 2: open flushed
        ser.reset_output_buffer()
        yield ser
    finally:
        ser.close()

def save_summary(save_dir, text, filename="summary_verbose.txt"):
    if save_dir is None:
        return None
    path = Path(save_dir) / filename
    path.write_text(text, encoding="utf-8")
    print(f"Saved {filename} to {path.as_posix()}")
    return path

def save_run_data(base_folder: str, data, filenames, sample_count, summary: str = ""):
    if LOAD_FROM_FILE:
        return None
    save_dir = get_save_dir(base_folder, sample_count)

    if len(data) != len(filenames):
        raise ValueError("Length of data and filenames must match.")
    for d, filename in zip(data, filenames):
        if d is not None:
            csv_path = save_dir / filename
            df = pd.DataFrame(d)
            df.to_csv(csv_path, index=False)
            print(f"Saved {len(d)} samples to {csv_path.as_posix()}")
        
    if summary:
        summary_path = save_dir / "summary.txt"
        summary_path.write_text(summary, encoding="utf-8")
        print(f"Saved summary to {summary_path.as_posix()}")
        
    return save_dir

def read_data_realtime(samples_total):
    # Fix 1: `with serial_port()` guarantees the port closes on any exit path.
    # Fix 2: opens flushed + drains before the loop.  Fix 3: try_read_realtime resyncs.
    chunks = []
    sample_count = 0
    try:
        with serial_port() as ser:
            drain_to_latest(ser)
            last_data = time.time()
            while samples_total is None or sample_count < samples_total:
                data = try_read_realtime(ser)
                if data is None:
                    if time.time() - last_data > 5:
                        print("No data for 5 s -- stopping capture.")
                        break
                    continue
                if np.mean(np.diff(data['timestamp'])) <= 0:
                    print("Timestamps did not change")
                    continue
                last_data = time.time()   # reset idle clock on a good packet
                chunks.append(data)
                sample_count += SAMPLES_PER_PACKET
    except KeyboardInterrupt:
        print("Stopping...")
    if not chunks:
        print("No data captured.")
        return np.empty(0, dtype=payload_dtype)   # always an ndarray -> no Optional warnings
    return np.concatenate(chunks)




# --- added when splitting the notebook ---
# In the notebook, Fs and CARRIER_FREQ were set by the calibration cell and shared through
# the kernel. As separate scripts, 01_calibrate.py saves them to calibration.json and every
# other script reads them here.
CALIBRATION_FILE = Path(__file__).resolve().parent / "calibration.json"
Fs = 1500.0          # same initial guess the calibration cell starts from
CARRIER_FREQ = 0.0

if CALIBRATION_FILE.exists():
    _cal = json.loads(CALIBRATION_FILE.read_text())
    Fs = float(_cal["Fs"])
    CARRIER_FREQ = float(_cal["CARRIER_FREQ"])
else:
    print(f"WARNING: {CALIBRATION_FILE.name} not found. Run 01_calibrate.py first. "
          f"Using Fs = {Fs} Hz and CARRIER_FREQ = {CARRIER_FREQ} Hz.")
