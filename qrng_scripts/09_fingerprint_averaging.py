r"""
## Fingerprint Averaging and Stability

**Purpose:** Repeatedly averages each sensor's complex carrier amplitude, turns every read into pairwise comparison bits, and reports how ratios and bitstrings change over time.

**Usage:** Set `AVERAGE_SAMPLES` to the samples coherently averaged in one read and `NUM_FINGERPRINTS` to the number of reads. The cell prints the required capture duration before acquisition. Larger averages reduce zero-mean noise but make each read slower. Keep `VERBOSE = False` for large runs.

For four sensors the bit order is A1>A2, A1>A3, A1>A4, A2>A3, A2>A4, A3>A4. The report shows per-sensor drift, pairwise-ratio drift, the modal bitstring, and how many distinct strings occurred. One distinct string means repeatability within this record; it does **not** prove uniqueness, entropy, or environmental robustness.

Live mode saves raw data, `summary.txt`, and figures under `data/averaging_data/...`. Saved-file mode uses the first file and reports if it contains fewer complete reads than requested.

Source: notebook cell 21 (code copied unchanged unless noted below).
"""

from qrng_common import *
from qrng_viz import *

### CONFIG

AVERAGE_SAMPLES  = 1    # samples averaged per fingerprint read -> a_hat is the mean over this
NUM_FINGERPRINTS = 1000000       # how many consecutive reads to take, to see whether the bits drift over time
SETTLE_TAU       = 20.0      # filter time constants to drop at the start (7*tau ~= -60 dB ring-down)
DEMOD_LP_HZ      = 20       # Hz, low-pass cutoff for the demodulated envelope
FOLDER           = "averaging_data"
VERBOSE          = False
PLOT_SUMMARY     = True     # render the summary figure (saved next to summary.txt)
PLOT_MAX_FINGERPRINTS = 10000000  # above this the per-read figure is unreadable -> skipped, with a notice

###

# capture length is DERIVED: settle the demod low-pass, then take the reads back to back.
if AVERAGE_SAMPLES < 1:
    raise ValueError("AVERAGE_SAMPLES must be at least 1")
TRANSIENT = int(np.ceil(SETTLE_TAU / (np.pi * DEMOD_LP_HZ) * Fs))    # filter roll-off, dropped below
NUM_SAMPLES_TO_SAVE = int(np.ceil((TRANSIENT + NUM_FINGERPRINTS * AVERAGE_SAMPLES)
                                  / SAMPLES_PER_PACKET) * SAMPLES_PER_PACKET)   # whole packets

print(f"{NUM_FINGERPRINTS} fingerprints x {AVERAGE_SAMPLES} samples ({AVERAGE_SAMPLES/Fs:.3f} s each) "
      f"+ {TRANSIENT} settling ({TRANSIENT/Fs:.3f} s)"
      f"  ->  capture {NUM_SAMPLES_TO_SAVE} samples ({NUM_SAMPLES_TO_SAVE/Fs:.1f} s)")

sos_lps = {s: signal.butter(4, DEMOD_LP_HZ, btype='low', fs=Fs, output='sos') for s in SENSORS}
zi_I   = {s: np.zeros((len(sos_lps[s]), 2)) for s in SENSORS}  # type: ignore
zi_Q   = {s: np.zeros((len(sos_lps[s]), 2)) for s in SENSORS}  # type: ignore


def filter_pipeline(data):
    # Returns the COMPLEX baseband Z = I + jQ.
    global sample_count
    N = SAMPLES_PER_PACKET
    t = np.arange(sample_count, sample_count + N) / Fs
    lo_cos, lo_sin = np.cos(2*np.pi*CARRIER_FREQ*t), -np.sin(2*np.pi*CARRIER_FREQ*t)
    sample_count += N          # pipeline owns LO phase for both live data and replay

    Z = {}
    for s in SENSORS:
        v_mv = adc_to_mv(data[s])
        I, zi_I[s] = signal.sosfilt(sos_lps[s], lo_cos * v_mv, zi=zi_I[s])
        Q, zi_Q[s] = signal.sosfilt(sos_lps[s], lo_sin * v_mv, zi=zi_Q[s])
        Z[s] = I + 1j * Q
    return Z


sample_count = 0
timeouts = 0
final_data = {s: [] for s in SENSORS}
raw_data = []

if LOAD_FROM_FILE:
    data = read_from_file_instant(LOAD_FILE_PATHS[0])[:NUM_SAMPLES_TO_SAVE]
    if data is not None:
        for i in range(len(data) // SAMPLES_PER_PACKET):
            z = filter_pipeline(data[i * SAMPLES_PER_PACKET:(i + 1) * SAMPLES_PER_PACKET])
            for s in SENSORS:
                final_data[s].append(z[s])
else:
    with serial_port() as ser:
        drain_to_latest(ser)
        while sample_count < NUM_SAMPLES_TO_SAVE:
            data = try_read_realtime(ser)
            if data is None:
                timeouts += 1
                if timeouts > 20:
                    print("serial timed out repeatedly - stopping early")
                    break
                continue
            timeouts = 0
            raw_data.append(data)

            z = filter_pipeline(data)  
            for s in SENSORS:     
                final_data[s].append(z[s])


if not final_data[SENSORS[0]]:
    raise RuntimeError("No complete packets were available for averaging.")
z_all = {s: np.concatenate(final_data[s])[TRANSIENT:] for s in SENSORS}
n     = len(z_all[SENSORS[0]])
t_all = np.arange(TRANSIENT, TRANSIENT + n) / Fs

# de-rotate
ref = max(SENSORS, key=lambda s: np.abs(z_all[s]).mean())
phi = np.unwrap(np.angle(z_all[ref]))
slope, offset = np.polyfit(t_all, phi, 1)
freq_err = slope / (2 * np.pi)                                  # Hz
sig = {s: z_all[s] * np.exp(-1j * (slope * t_all + offset)) for s in SENSORS}

n_fp = min(NUM_FINGERPRINTS, n // AVERAGE_SAMPLES)
if n_fp < 1:
    raise RuntimeError(f"only {n} samples after settling - not enough for one "
                       f"{AVERAGE_SAMPLES}-sample read")
if n_fp < NUM_FINGERPRINTS:
    print(f"Short record: {n_fp} complete fingerprint(s) of the {NUM_FINGERPRINTS} requested "
          f"({n} settled samples available).")

a_hat, bits = {}, {}
for g in range(n_fp):
    sl = slice(g * AVERAGE_SAMPLES, (g + 1) * AVERAGE_SAMPLES)
    a_hat[g] = {s: 2 * np.abs(np.mean(sig[s][sl])) for s in SENSORS}
    bits[g] = [1 if a_hat[g][SENSORS[i]] > a_hat[g][SENSORS[j]] else 0
               for i in range(len(SENSORS)) for j in range(i + 1, len(SENSORS))]

# --- report ------------------------------------------------------------------------
result  = f"LO frequency error   df = {freq_err:+.4f} Hz  (large |df| -> re-run Setup / refine CARRIER_FREQ)\n"
result += (f"{n_fp} fingerprint read(s) x {AVERAGE_SAMPLES} samples averaged ({AVERAGE_SAMPLES/Fs:.3f} s each), "
           f"{n} samples after settling ({TRANSIENT} dropped, {TRANSIENT/Fs:.3f} s)\n")
if n_fp < NUM_FINGERPRINTS:
    result += f"NOTE: {NUM_FINGERPRINTS} reads were requested; the record only held {n_fp}.\n"
result += "\n"

if VERBOSE:
    for g in range(n_fp):
        means = "  ".join(f"{s}={a_hat[g][s]:7.3f}" for s in SENSORS)
        result += f"Read {g+1:2d}: {means}   bits={bits[g]}\n"

result += "\nPer-sensor drift across reads (mV) - the floor the surface fingerprint must beat:\n"
for s in SENSORS:
    v = np.array([a_hat[g][s] for g in range(n_fp)])
    result += f"  {s:9s} mean={v.mean():7.3f}  std={v.std():7.4f}  CV={100*v.std()/abs(v.mean()):5.2f}%\n"

result += "\nPairwise ratio drift (common mode cancels in the ratio):\n"
for i in range(len(SENSORS)):
    for j in range(i + 1, len(SENSORS)):
        r = np.array([a_hat[g][SENSORS[i]] / a_hat[g][SENSORS[j]] for g in range(n_fp)])
        result += f"  {SENSORS[i]:9s}/{SENSORS[j]:9s} mean={r.mean():7.4f}  std={r.std():7.4f}  CV={100*r.std()/abs(r.mean()):5.2f}%\n"

allbits = np.array([bits[g] for g in range(n_fp)])
uniq, counts = np.unique(allbits, axis=0, return_counts=True)
modal_bitstring = [int(x) for x in uniq[counts.argmax()]]
result += f"\nBit stability: {counts.max()}/{n_fp} reads gave the bitstring {' '.join(map(str, modal_bitstring))}\n"
result += f"{len(uniq)} distinct bitstring(s) seen (1 = perfectly reproducible)\n"

print(result)

raw_out = np.concatenate(raw_data) if raw_data else None
save_dir = save_run_data(FOLDER, data=[raw_out], filenames=["data.csv"], sample_count=NUM_SAMPLES_TO_SAVE, summary=result)

# summary figure (bit stability, per-read drift in ppm, per-sensor vs pairwise CV)
if PLOT_SUMMARY:
    if n_fp > PLOT_MAX_FINGERPRINTS:
        print(f"Summary figure skipped: {n_fp} reads > PLOT_MAX_FINGERPRINTS ({PLOT_MAX_FINGERPRINTS}).")
    else:
        plot_averaging_summary(a_hat, bits, freq_err, AVERAGE_SAMPLES, save_dir=save_dir, trend_points=100)
        plot_averaging_table_summary(a_hat, bits, freq_err, AVERAGE_SAMPLES, save_dir=save_dir)
