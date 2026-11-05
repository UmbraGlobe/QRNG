r"""
## Online Common-Mode Removal (RLS)

**Purpose:** Uses one sensor as an adaptive reference and subtracts the part of each other sensor that follows it. The outputs are experimental common-mode-rejected residuals.

**Usage:** Choose `REF` from `SENSORS`, set `NUM_SAMPLES_TO_SAVE`, and run after calibration. The reference channel is not included in the output. `BANDWIDTH` controls the settling guard and should remain consistent with the ±2 Hz carrier filter used here.

Live mode saves one processed CSV per non-reference sensor under `data/rls_filtered_data/...`. These residual files do not contain raw timestamps and ADC columns, so do not use them as normal `LOAD_FILE_PATHS` input. Saved-file mode runs the filter without creating a new run folder.

Source: notebook cell 19 (code copied unchanged unless noted below).
"""

from qrng_common import *

### CONFIG

NUM_SAMPLES_TO_SAVE = 20000 
FOLDER = "rls_filtered_data"
REF = SENSORS[0]
BANDWIDTH = 4
SETTLE_TAU = 7.0 

###

USABLE_SENSORS = [s for s in SENSORS if s != REF]
TRANSIENT = int(np.ceil(SETTLE_TAU / (np.pi * BANDWIDTH) * Fs)) 
TRANSIENT = int(np.ceil(TRANSIENT/1000) * 1000)  # round up to nearest 1000 samples

sos_bps = {s: signal.butter(4, [CARRIER_FREQ - 2, CARRIER_FREQ + 2], btype='bp', fs=Fs, output='sos') for s in SENSORS}
sos_lps = {s: signal.butter(4, 20, btype='low', fs=Fs, output='sos') for s in SENSORS}
zi_bps = {s: np.zeros((len(sos_bps[s]), 2)) for s in SENSORS} # type: ignore
zi_I  = {s: np.zeros((len(sos_lps[s]), 2)) for s in SENSORS} # type: ignore
zi_Q  = {s: np.zeros((len(sos_lps[s]), 2)) for s in SENSORS} # type: ignore
sample_count = 0

NTAPS = 1
filters  = {s: pa.filters.FilterRLS(n=NTAPS, mu=0.99999) for s in USABLE_SENSORS} # type: ignore
ref_hist = np.zeros(NTAPS)       # persistent rolling window of A[REF] across packets

def filter_pipeline(data):
    global sample_count, ref_hist
    N = SAMPLES_PER_PACKET
    t = np.arange(sample_count, sample_count + N) / Fs      # continuous LO phase
    lo_cos, lo_sin = np.cos(2*np.pi*CARRIER_FREQ*t), -np.sin(2*np.pi*CARRIER_FREQ*t)
    sample_count += N
 
    A = {}
    for s in SENSORS:
        v_mv = adc_to_mv(data[s])
        v_bp, zi_bps[s] = signal.sosfilt(sos_bps[s], v_mv, zi=zi_bps[s])
        I, zi_I[s] = signal.sosfilt(sos_lps[s], lo_cos*v_bp, zi=zi_I[s])
        Q, zi_Q[s] = signal.sosfilt(sos_lps[s], lo_sin*v_bp, zi=zi_Q[s])
        A[s] = 2.0 * np.sqrt(I**2 + Q**2)

    errors = {s: np.empty(N) for s in USABLE_SENSORS}
    ref = A[REF]
    for i in range(N):
        x = ref[i:i+1]                        # length-1 input window
        for s in USABLE_SENSORS:
            y = filters[s].predict(x)         # estimated common mode = w·A_REF
            filters[s].adapt(A[s][i], x)      # desired = sensor sample
            errors[s][i] = A[s][i] - y        # cleaned = sensor − common mode
    return errors

final_data = {s: [] for s in USABLE_SENSORS}
reader = None
timeout = 0
sample_count = 0
if LOAD_FROM_FILE and reader is None:
    reader = read_from_file_realtime(LOAD_FILE_PATHS[0], pace=False)
with serial_port() as ser:                      # Fix 1: guaranteed close (no-op when LOAD_FROM_FILE)
    if not LOAD_FROM_FILE:
        drain_to_latest(ser)
    while NUM_SAMPLES_TO_SAVE is None or sample_count < NUM_SAMPLES_TO_SAVE:
        data = reader() if LOAD_FROM_FILE and reader is not None else try_read_realtime(ser)
        if data is None:
            if LOAD_FROM_FILE and reader is not None and reader.exhausted:
                break
            continue
        d = filter_pipeline(data)
        for s in USABLE_SENSORS:
            final_data[s].append(d[s])


if sample_count <= TRANSIENT:
    raise RuntimeError(f"Only {sample_count} samples were processed; more than {TRANSIENT} are required.")
sig = {s: np.concatenate(final_data[s])[TRANSIENT:] for s in USABLE_SENSORS}

save_run_data(FOLDER, data=[sig[s] for s in USABLE_SENSORS], filenames=[f"data_{s}.csv" for s in USABLE_SENSORS], sample_count=sample_count, summary=f"Sampling Rate = {Fs:.2f} Hz\nCarrier Frequency = {CARRIER_FREQ:.6f} Hz\n")
