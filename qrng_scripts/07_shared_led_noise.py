r"""
## Shared LED Fluctuation and Independent Sensor Noise

**Purpose:** Compares all sensor envelopes to separate fluctuation shared by the LED from noise local to each photodiode channel.

The model is $A_i(t)=a_i+a_i m(t)+n_i(t)$, where $a_i$ is mean carrier amplitude, $m(t)$ is common fractional LED fluctuation, and $n_i(t)$ is independent channel noise. Cross-covariance estimates the shared term because independent noise should not correlate between channels:

$$\operatorname{Var}(m)=\frac{\operatorname{Cov}(A_i,A_j)}{a_i a_j}, \qquad \operatorname{Var}(n_i)=\operatorname{Var}(A_i)-a_i^2\operatorname{Var}(m).$$

**Usage:** Choose `NUM_TRIALS` and `NUM_SAMPLES_TO_SAVE` per trial. Saved-file mode accepts one long file or exactly one file per trial. Keep `CORR_GATE = 0.5` unless intentionally changing which pairs count as coupled. For the standard ±2 Hz working band, use `BANDWIDTH = 4`; treat the current experimental `0.1` setting cautiously because extremely narrow filters can create false common mode. Use about 60 seconds of usable data per trial when quoting `RMS(m)`.

**Interpretation:** `corr_A` near 1 means channels move together; `RMS(m) = 0.01` means 1% shared fractional fluctuation; `Var(n_i)` is independent noise power in `mV²`. `NaN Var(m)` means no pair passed the correlation gate. Correlation of `1.0000` together with zero independent noise is a warning for filter transient or another artifact, not proof of ideal behavior.

Live mode saves raw trials, `summary.txt`, and figures under `data/cross_sensor_variance_data/...`. The current untextured A2V-16 array provides instrument measurements, not PUF uniqueness results.

Source: notebook cell 17 (code copied unchanged unless noted below).
"""

from qrng_common import *
from qrng_viz import *

### CONFIG

NUM_SAMPLES_TO_SAVE = 1000000
NUM_TRIALS = 1
CORR_GATE = 0.5      # only average pairs that actually share the common mode
FOLDER = "cross_sensor_variance_data"
BANDWIDTH = 0.1       # Hz, width of the bandpass filter around the carrier frequency
DEMOD_LP_HZ = 20    # Hz, demod low-pass cutoff
SETTLE_TAU = 7.0    
PLOT_SUMMARY = True  # render the per-trial summary figure (saved next to summary.txt)

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

# filters:
# sos_bp: isolate the carrier (CARRIER_FREQ +/- BW/2 Hz)
# sos_lp: extract the demodulated envelope after mixing the carrier down to baseband
# Both zero-phase (filtfilt)
sos_bp = signal.butter(4, [CARRIER_FREQ - BANDWIDTH/2, CARRIER_FREQ + BANDWIDTH/2], btype='bp', fs=Fs, output='sos')
sos_lp = signal.butter(4, DEMOD_LP_HZ, btype='low', fs=Fs, output='sos')
offdiag = ~np.eye(len(SENSORS), dtype=bool)   # mask that selects off-diagonal (cross) entries

TRANSIENT_BP = int(np.ceil(SETTLE_TAU / (np.pi * BANDWIDTH) * Fs))
TRANSIENT_LP = int(np.ceil(SETTLE_TAU / (np.pi * DEMOD_LP_HZ) * Fs))
TRANSIENT_BP = int(np.ceil(TRANSIENT_BP / 1000) * 1000)  # round up to nearest 1000 samples
TRANSIENT_LP = int(np.ceil(TRANSIENT_LP / 1000) * 1000)  # round up to nearest 1000 samples


result = "A_i is the demodulated message of sensor i (LED FLUCTUATION + CIRCUIT NOISE + INTERFERENCE...)\n"
result += "a_i is the mean carrier amplitude of sensor i\n"
result += "Cov(A_i, A_j) = a_i * a_j * Var(m)  ->  Var(m) = Cov(A_i, A_j) / (a_i * a_j)\n"
result += "Var(n_i) = Var(A_i) - a_i^2 * Var(m)  (independent per-diode noise)\n"
result += f"\nAcross {len(SENSORS)} sensors over {NUM_TRIALS} independent windows, Fs = {Fs:.2f} Hz, Carrier Frequency = {CARRIER_FREQ:.2f} Hz\n"

var_m_trials = []
trial_stats = {}     # ADDED: per-trial numbers kept for the summary figure

for trial in range(NUM_TRIALS):
    # raw ADC -> mV -> bandpass around the carrier, then drop the filter transient
    bandpassed = {}
    for s in SENSORS:
        v_mv = adc_to_mv(data[trial][s])
        v_bp = signal.sosfiltfilt(sos_bp, v_mv)
        bandpassed[s] = v_bp[TRANSIENT_BP:-TRANSIENT_BP]

    n = min(len(bandpassed[s]) for s in SENSORS)   # common length across channels
    t = np.arange(n) / Fs

    # 1. Complex I/Q Downconversion
    lo_cos = np.cos(2 * np.pi * CARRIER_FREQ * t)
    lo_sin = -np.sin(2 * np.pi * CARRIER_FREQ * t)

    Z = {}
    for s in SENSORS:
        I = signal.sosfiltfilt(sos_lp, lo_cos * bandpassed[s][:n])
        Q = signal.sosfiltfilt(sos_lp, lo_sin * bandpassed[s][:n])
        Z[s] = I + 1j * Q

    # 2. De-rotate using highest-amplitude channel as phase reference
    ref = max(SENSORS, key=lambda s: np.abs(Z[s]).mean())
    phi = np.unwrap(np.angle(Z[ref]))
    slope, offset = np.polyfit(t, phi, 1)
    freq_err = slope / (2 * np.pi)  # Carrier frequency error in Hz

    # Apply phase correction vector across all channels
    sig = {s: Z[s] * np.exp(-1j * (slope * t + offset)) for s in SENSORS}

    # 3. Recover each channel's de-rotated envelope A_i(t) and mean amplitude a_i
    A = {}
    a = {}
    for s in SENSORS:
        envelope = 2.0 * np.abs(sig[s])
        A[s] = envelope[TRANSIENT_LP:-TRANSIENT_LP]  # drop low-pass edge transient
        a[s] = A[s].mean()

    X = np.vstack([A[s] for s in SENSORS])
    a_vec  = np.array([a[s] for s in SENSORS])
    cov_A  = np.cov(X) # mV^2, Var(A_i) on the diagonal
    corr_A = np.corrcoef(X)  # dimensionless, in [-1, 1]
    norm_cov = cov_A / np.outer(a_vec, a_vec) # Cov(A_i,A_j) / (a_i * a_j)

    # mask out entries that aren't correlated; var(m) = mean of what survives (nan if none)
    gated = np.where(offdiag & (corr_A > CORR_GATE), norm_cov, np.nan)
    var_m = np.nanmean(gated)
    n_pairs = int(np.isfinite(gated).sum() // 2)
    is_coupled = np.isfinite(gated).sum(axis=1) > 0
    var_m_trials.append(var_m)

    var_A = np.diag(cov_A)  # total envelope variance per channel
    var_common = a_vec ** 2 * var_m # portion explained by the shared LED
    var_n = var_A - var_common # remainder = independent, uncorrelated noise

    trial_stats[PERIODS[trial]] = dict(cov_A=cov_A, corr_A=corr_A, norm_cov=norm_cov,
                                      a_vec=a_vec, var_A=var_A, var_common=var_common,
                                      var_n=var_n, var_m=var_m, n_pairs=n_pairs,
                                      freq_err=freq_err)

    result += f"\n{'=' * 64}\nTRIAL {PERIODS[trial]}\n{'=' * 64}\n"
    result += f"LO frequency error df = {freq_err:+.4f} Hz (ref channel: {ref})\n\n"
    result += fmt_matrix(cov_A, "[1] COVARIANCE (mV^2) - diagonal = total fluctuations of each channel, off-diagonal = fluctuation shared between two channels")
    result += fmt_matrix(corr_A, "[2] CORRELATION (-1 to 1)- near 1 = channels track togther (shared LED), 0 = independent")
    result += fmt_matrix(norm_cov, "[3] NORMALIZED COVARIANCE - off-diagonals = var(m) estimate per pair; all equal => one shared common mode fluctuation", fmt="{:>12.4e}")
    result += f"\nvar(m) = {var_m:.4e}   RMS(m) = sqrt(var(m)) = {np.sqrt(var_m):.4f}"
    result += f"   (mean over {n_pairs} coupled pair(s), corr > {CORR_GATE})\n"

    # per-channel noise split: Var(A_i) = a_i^2*Var(m) [common] + Var(n_i) [independent]
    result += f"\n{'channel':>10}{'a_i (mV)':>10}{'Var(A_i)':>12}{'a^2*Var(m)':>12}{'Var(n_i)':>12}{'% indep':>9}{'coupled':>9}\n"
    clamped_any = False
    for i, s in enumerate(SENSORS):
        vn = var_n[i]
        if vn < 0: # unphysical: a_i^2*Var(m) > Var(A_i)
            vn = 0.0  # catastrophic cancellation -> report as ~0
            clamped_any = True
        pct_indep = 100 * vn / var_A[i]
        flag = "yes" if is_coupled[i] else "NO"
        result += (f"{s:>10}{a_vec[i]:>10.3f}{var_A[i]:>12.4f}" f"{var_common[i]:>12.4f}{vn:>12.4f}{pct_indep:>8.1f}%{flag:>9}\n")
    if clamped_any:
        result += "  Var(n_i)=0 (clamped): a_i^2*Var(m) >= Var(A_i), i.e. ~100% common mode; independent noise is below this method's resolution.\n"

result += f"\nE[RMS(m)] over {NUM_TRIALS} trial(s) = {np.sqrt(np.nanmean(var_m_trials)):.4f}  (fractional LED intensity fluctuation)\n"

result += f"\n{'=' * 64}\nINTERPRETATION\n{'=' * 64}\n"
result += (
    "Here are the variables in the model V(t) = A_i(t)cos(2πf_c t + θ),  A_i(t) = a_i + a_i·m(t) + n_i(t):\n\n"
    "Measured / raw\n"
    "    V_i(t)      — the raw load voltage of diode i (the carrier sine you actually sample), mV.\n"
    "    f_c         — carrier frequency.\n\n"
    "Per-diode envelope (after I/Q demod & de-rotation)\n"
    "    A_i(t)      — demodulated envelope of diode i: the slowly-varying carrier amplitude, mV. "
    "This is what all the covariance math runs on.\n"
    "    a_i         — mean carrier amplitude of diode i (the constant part of A_i), mV. "
    "Represents the average optical power reaching that diode.\n\n"
    "The two things we're separating\n"
    "    m(t)        — common mode fractional intensity fluctuation, dimensionless, E[m]=0. "
    "Shared by all diodes (common mode).\n"
    "    n_i(t)      — independent per-diode in-band noise (shot/thermal + local pickup), mV. "
    "Uncorrelated across diodes. Currently pinned at the electronics floor.\n\n"
    "Computed statistics\n"
    "    Var(m)      — variance of m, dimensionless. Estimated as Cov(A_i,A_j)/(a_i·a_j).\n"
    "    RMS(m)      — √Var(m), the fractional fluctuation (0.01 = 1%). E[RMS(m)] averages it over trials.\n"
    "    Var(A_i)    — total variance of diode i's envelope, mV² (diagonal of cov_A).\n"
    "    a_i²·Var(m) — the part of Var(A_i) explained by the shared LED (common mode), mV².\n"
    "    Var(n_i)    — remainder Var(A_i) − a_i²·Var(m) = independent noise, mV².\n"
    "    cov_A       — envelope covariance matrix, mV² (Var(A_i) on diagonal, shared fluctuation off-diagonal).\n"
    "    corr_A      — envelope correlation matrix, −1…1 (near 1 = tracking the same LED).\n"
    "    norm_cov    — cov_A / (a_i·a_j); off-diagonals are the per-pair Var(m) estimate.\n"
)

print(result)

data_to_save = [data[i] for i in range(NUM_TRIALS)]
file_names = [f"data_{PERIODS[i]}.csv" for i in range(NUM_TRIALS)]
save_dir = save_run_data(FOLDER, data=data_to_save, filenames=file_names, sample_count=NUM_SAMPLES_TO_SAVE, summary=result)

if PLOT_SUMMARY:
    for lbl, d in trial_stats.items():
        plot_cross_sensor_summary(**d, label=lbl, save_dir=save_dir)
