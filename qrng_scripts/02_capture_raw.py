r"""
## Capture and Save Raw Trials

**Purpose:** Captures one uninterrupted live stream and splits it into equal raw CSV files for replay or later comparison.

**Usage:** This section is live-capture only. Set `NUM_SAMPLES_TO_SAVE` to the samples in each trial and `NUM_TRIALS` to the number of files required. Total capture size is their product. Run the cell and check the printed saved-sample counts; an interrupted stream can leave the last trial short.

Files are saved as `data_0.csv`, `data_1.csv`, and so on under `data/raw_data/...`. If `LOAD_FROM_FILE = True`, capture is deliberately skipped.

Source: notebook cell 7 (code copied unchanged unless noted below).
"""

from qrng_common import *

NUM_SAMPLES_TO_SAVE = 3000
NUM_TRIALS = 10
data = {}
if LOAD_FROM_FILE:
    print("LOAD_FROM_FILE is enabled; raw capture is skipped.")
else:
    d = read_data_realtime(NUM_SAMPLES_TO_SAVE * NUM_TRIALS)
    for i in range(NUM_TRIALS):
        data[i] = d[i * NUM_SAMPLES_TO_SAVE:(i + 1) * NUM_SAMPLES_TO_SAVE]
    save_run_data("raw_data", data=[data[i] for i in range(NUM_TRIALS)], sample_count=NUM_SAMPLES_TO_SAVE, filenames=[f"data_{i}.csv" for i in range(NUM_TRIALS)])