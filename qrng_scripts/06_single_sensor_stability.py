r"""
## Single-Sensor Stability Across Time Windows

**Purpose:** Tests whether one sensor's mean demodulated carrier amplitude and total fluctuation remain stable across consecutive time windows.

**Usage:** Choose `SENSOR`, `NUM_TRIALS`, and `NUM_SAMPLES_TO_SAVE` per trial. Saved-file mode accepts either one continuous file containing all trials or exactly one file per trial. Every trial must contain the requested samples. Keep the standard filter settings unless deliberately testing bandwidth; filter-edge samples are discarded automatically with `SETTLE_TAU`.

The report compares mean amplitude `a`, envelope variance `Var(A)`, standard deviation, and fractional fluctuation `std/a`. The coefficient of variation (CV) shows how much `Var(A)` changes between windows; smaller CV means a more stationary noise level. This measures within-run stability, not power-cycle or day-to-day drift.

Live mode saves each trial, `summary.txt`, and a summary plot under `data/temporal_variance_data/...`.

Source: notebook cell 15 (code copied unchanged unless noted below).
"""

from qrng_common import *
from qrng_viz import *


### CONFIG

NUM_SAMPLES_TO_SAVE = 50000   # set to None to run until interrupted
FOLDER = "temporal_variance_data"
SENSOR = SENSORS[0]
NUM_TRIALS = 5
DEMOD_LP_HZ = 20  # Hz, low-pass filter cutoff for demodulated signal
TRANSIENT = 200  # samples to drop to ensure filter settles
SETTLE_TAU = 7.0  # filter time constants to drop at the start (7*tau ~= -60 dB ring-down)
BANDWIDTH = 4  # Hz, width of the bandpass filter around the carrier frequency
PLOT_SUMMARY = True  # render the summary figure (saved next to summary.txt)

###

PERIODS = [f"t{i+1}" for i in range(NUM_TRIALS)]
if LOAD_FROM_FILE and len(LOAD_FILE_PATHS) not in (1, NUM_TRIALS):
    raise ValueError(f"Provide one continuous file or {NUM_TRIALS} per-trial files.")

data = {}   
if LOAD_FROM_FILE:
    if len(LOAD_FILE_PATHS) == 1:
        d = read_from_file_instant(LOAD_FILE_PATHS[0])[:NUM_SAMPLES_TO_SAVE * NUM_TRIALS]
        data = {i: d[i * NUM_SAMPLES_TO_SAVE:(i + 1) * NUM_SAMPLES_TO_SAVE]
                for i in range(NUM_TRIALS)}
    else:
        data = {i: read_from_file_instant(LOAD_FILE_PATHS[i])[:NUM_SAMPLES_TO_SAVE]
                for i in range(NUM_TRIALS)}
else:
    d = read_data_realtime(NUM_SAMPLES_TO_SAVE * NUM_TRIALS)
    for i in range(NUM_TRIALS):
        data[i] = d[i * NUM_SAMPLES_TO_SAVE:(i + 1) * NUM_SAMPLES_TO_SAVE]

short = {PERIODS[i]: len(data[i]) for i in range(NUM_TRIALS)
         if len(data[i]) < NUM_SAMPLES_TO_SAVE}
if short:
    raise RuntimeError(f"Incomplete trial data: {short}; expected {NUM_SAMPLES_TO_SAVE} each.")

TRANSIENT_BP = int(np.ceil(SETTLE_TAU / (np.pi * BANDWIDTH) * Fs))
TRANSIENT_LP = int(np.ceil(SETTLE_TAU / (np.pi * DEMOD_LP_HZ) * Fs))
TRANSIENT_BP = int(np.ceil(TRANSIENT_BP / 1000) * 1000)  # round up to nearest 1000 samples
TRANSIENT_LP = int(np.ceil(TRANSIENT_LP / 1000) * 1000)  # round up to nearest 1000 samples


data_1sensor = {i: data[i][SENSOR] for i in range(NUM_TRIALS)}
        
sos_bp = signal.butter(4, [CARRIER_FREQ - BANDWIDTH/2, CARRIER_FREQ + BANDWIDTH/2], btype='bp', fs=Fs, output='sos')
sos_lp = signal.butter(4, DEMOD_LP_HZ, btype='low', fs=Fs, output='sos')

filtered_data = {i: signal.sosfiltfilt(sos_bp, data_1sensor[i]) for i in range(NUM_TRIALS)}
kept_filtered_data = {i: filtered_data[i][TRANSIENT_BP:-TRANSIENT_BP] for i in range(NUM_TRIALS)}


n = min(len(kept_filtered_data[i]) for i in range(NUM_TRIALS))
t = np.arange(n) / Fs 

A = {}
a = {}
for i in range(NUM_TRIALS):
    x = kept_filtered_data[i][:n]
    I = np.cos(2 * np.pi * CARRIER_FREQ * t) * x
    Q = -np.sin(2 * np.pi * CARRIER_FREQ * t) * x
    I_f = signal.sosfiltfilt(sos_lp, I)
    Q_f = signal.sosfiltfilt(sos_lp, Q)
    A[i] = (2 * np.sqrt(I_f**2 + Q_f**2))[TRANSIENT_LP:-TRANSIENT_LP]
    a[i] = np.mean(A[i])

a_vec = np.array([a[i] for i in range(NUM_TRIALS)])

var_A = np.array([np.var(A[i]) for i in range(NUM_TRIALS)])  # mV^2, total fluctuation per window
std_A = np.sqrt(var_A) # mV
frac  = std_A / a_vec # std/mean: fractional fluctuation (m + n combined)

result = f"A is the message encoded in the carrier wave (the total noise signal, LED FLUCTUATION + CIRCUIT NOISE + INTERFERENCE...)\n" 
result += f"a is the mean of the demodulated message\n"
result += f"Var(A) is the total variance of A\n"
result += f"std(A) is the standard deviation of A\n"
result += f"std/a is the fractional fluctuation (std/mean) for each time window.\n"

result += f"\nSame sensor ({SENSOR}) over {NUM_TRIALS} sequential time windows, Fs = {Fs:.2f} Hz\n"
result += f"{'window':>8}{'a (mV)':>12}{'Var(A) mV^2':>14}{'std(A) mV':>12}{'std/a':>12}\n"
for i in range(NUM_TRIALS):
    result += f"{PERIODS[i]:>8}{a_vec[i]:>12.4f}{var_A[i]:>14.4f}{std_A[i]:>12.4f}{frac[i]:>12.4e}\n"

# across-window summary: the spread of Var(A) shows how stationary the total
# fluctuation is. CV = std/mean of Var(A) across windows; small => stationary.
mean_var = var_A.mean()
std_var  = var_A.std(ddof=1)
cv = std_var / mean_var

print(result)
result += f"\nVar(A) across windows: E[Var(A)] = {mean_var:.4f} mV^2, std = {std_var:.4f} mV^2, spread (CV) = {100*cv:.1f}%"

data_to_save = [data[i] for i in range(NUM_TRIALS)]
file_names = [f"data_{PERIODS[i]}.csv" for i in range(NUM_TRIALS)]
save_dir = save_run_data(FOLDER, data=data_to_save, filenames=file_names, sample_count=NUM_SAMPLES_TO_SAVE, summary=result)   # CHANGED: keep the run folder


if PLOT_SUMMARY:
    plot_summary(f'Window variance — {SENSOR}', [
        dict(kind='stat', span=2, items=[
            ('E[Var(A)]', f'{mean_var:.4f} mV²', f'across {NUM_TRIALS} windows'),
            ('std of Var(A)', f'{std_var:.4f} mV²', 'window-to-window'),
            ('spread (CV)', f'{100*cv:.1f} %', 'small => stationary',
             'good' if cv < 0.10 else 'warning'),
            ('mean a', f'{a_vec.mean():.1f} mV',
             f'spread {a_vec.min():.1f}–{a_vec.max():.1f} mV'),
            ('worst std/a', f'{frac.max():.2e}', 'fractional fluctuation (m + n)'),
        ]),
        dict(kind='bars', labels=PERIODS, values=a_vec, ylabel='a  (mV)', fmt='{:.2f}',
             title='Mean carrier amplitude per window'),
        dict(kind='bars', labels=PERIODS, values=var_A, ylabel='Var(A)  (mV²)', fmt='{:.4f}',
             ref=mean_var, ref_label='E[Var(A)]',
             title='Total envelope variance per window'),
        dict(kind='table', span=2, height=0.70 + 0.40 * (NUM_TRIALS + 1),
             title='Per-window statistics',
             columns=['window', 'a (mV)', 'Var(A) mV²', 'std(A) mV', 'std/a'],
             rows=[[PERIODS[i], f'{a_vec[i]:.4f}', f'{var_A[i]:.4f}',
                    f'{std_A[i]:.4f}', f'{frac[i]:.4e}'] for i in range(NUM_TRIALS)]),
    ], subtitle=f'Fs = {Fs:.0f} Hz   carrier = {CARRIER_FREQ:.1f} Hz   '
                f'{NUM_SAMPLES_TO_SAVE} samples/window',
       save_dir=save_dir, filename=f'window_variance_{SENSOR}.png')
