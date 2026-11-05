r"""
## Per-Channel Noise Floor and SNR

**Purpose:** Estimates each channel's independent broadband noise from carrier-free spectral bands, scales it to the working demodulation bandwidth, and compares it with carrier power.

**Usage:** Use a long, normally lit record and set `NUM_SAMPLES_TO_SAVE`; 500,000 samples take about 100 seconds at 5 kHz. `NFFT` controls frequency resolution. Keep the comb guards unless hardware timing changed, and keep `DEMOD_NOISE_BW` consistent with the total noise bandwidth used when comparing `Var(n_i)` elsewhere.

`noise uV rms` is standard deviation in the chosen band, `Var(n)` is the same noise power in `mV²`, and `SNR dB` compares carrier power with that noise. `consistent` means clean bands agree; `varies*` indicates frequency-dependent noise; `no clean band found` means the current gates could not isolate a trustworthy floor.

Live mode saves raw data, `summary.txt`, `summary_verbose.txt`, and a figure under `data/noise_floor_data/...`. Because this method avoids the carrier band, use the Dark Test to measure phase-locked pickup at the carrier.

Source: notebook cell 27 (code copied unchanged unless noted below).
"""

from qrng_common import *
from qrng_viz import *

### CONFIG

NUM_SAMPLES_TO_SAVE = 500000       # samples to load
NFFT                = 8192         # Welch segment length -> bin width = Fs/NFFT
BAND_HALFWIDTHS     = [2, 5, 10, 20]   # Hz, swept per band; a true floor is bandwidth-invariant (density)
DEMOD_NOISE_BW      = 4.0          # Hz, analysis band (±2 Hz bandpass) the floor is scaled to for Var(n)
CARRIER_HW          = 4.0          # Hz, half-width to integrate the carrier tone (signal power for the SNR)
GUARD_HZ            = 3.0          # keep bands this far from every known comb line
LOW_CUT_HZ         = 20.0         # ignore DC / drift / 1-f below this
PEAK_DB            = 6.0          # a band is contaminated if any bin exceeds its median by > this
DUTY               = 0.50         # PWM duty, to exclude (future) square-wave harmonic aliases too
HARMONIC_DBC       = -50.0        # only exclude carrier harmonics stronger than this vs the fundamental
FORCE_BANDS_HZ     = []           # extra band centers to check even if the auto-scan skips them
FOLDER             = "noise_floor_data"
PLOT_SUMMARY       = True         


###

nyq    = Fs / 2.0
bin_hz = Fs / NFFT

# one long record, same read path as the other analysis cells
data   = read_from_file_instant(LOAD_FILE_PATHS[0]) if LOAD_FROM_FILE else read_data_realtime(NUM_SAMPLES_TO_SAVE)
data   = data[:NUM_SAMPLES_TO_SAVE]
n_samp = len(data)
nper   = min(NFFT, n_samp)

# Welch PSD per sensor (mV^2/Hz)
psd = {}
for s in SENSORS:
    freqs, psd[s] = welch(adc_to_mv(data[s]), fs=Fs, nperseg=nper, noverlap=nper // 2,
                          window='hann', detrend='constant')

# ---- forbidden comb: every line a "carrier-free" band must avoid -------------
def fold(fr):                                 # alias a tone into [0, nyq]
    fr = abs(fr) % Fs
    return Fs - fr if fr > nyq else fr

def pwm_dbc(k):                               # k-th harmonic of a duty-DUTY pulse train, dBc vs fundamental
    ak, a1 = abs(np.sinc(k * DUTY)), abs(np.sinc(DUTY))
    return -np.inf if ak == 0 else 20 * np.log10(ak / a1)

forbidden = [float(CARRIER_FREQ)]                                       # carrier fundamental
for k in range(2, 200):                                                 # + (future PWM) harmonic aliases
    if pwm_dbc(k) >= HARMONIC_DBC:
        forbidden.append(float(fold(k * CARRIER_FREQ)))
forbidden += [float(500 * m) for m in range(1, int(nyq // 500) + 1)]    # 1 kHz MCU event -> 500/1000/... lines
forbidden += [float(50  * m) for m in range(1, int(nyq // 50)  + 1)]    # serial packet comb
forbidden  = np.array(sorted({round(x, 3) for x in forbidden if x > 0}))

FORBIDDEN_RANGES = [(400.0, 460.0)]           # wandering non-clock-locked interferer (hardware.md)
def in_forbidden_range(c):
    return any(rlo <= c <= rhi for rlo, rhi in FORBIDDEN_RANGES)

# ---- candidate bands = midpoints of the widest clean gaps in the comb --------
edges = np.array(sorted({LOW_CUT_HZ, nyq} | set(forbidden[(forbidden > LOW_CUT_HZ) & (forbidden < nyq)])))
cand = []
for a, b in zip(edges[:-1], edges[1:]):
    c = 0.5 * (a + b)
    hw_avail = 0.5 * (b - a) - GUARD_HZ
    if hw_avail >= min(BAND_HALFWIDTHS) and not in_forbidden_range(c):
        cand.append((round(c, 1), min(hw_avail, max(BAND_HALFWIDTHS))))
for x in FORCE_BANDS_HZ:
    cand.append((float(x), max(BAND_HALFWIDTHS)))
cand = sorted(set(cand))

# ---- measure the floor in each band, per sensor ------------------------------
def band_stats(pxx, center, hw):
    seg = pxx[(freqs >= center - hw) & (freqs <= center + hw)]
    if len(seg) < 4:
        return None
    return float(np.median(seg)), 10 * np.log10(seg.max() / np.median(seg))   # density mV^2/Hz, worst-bin dB

records = {s: [] for s in SENSORS}
for center, hw_band in cand:
    hw_top = min(hw_band, max(BAND_HALFWIDTHS))
    for s in SENSORS:
        dens_by_hw = []
        for hw in BAND_HALFWIDTHS:
            if hw > hw_band:
                continue
            st = band_stats(psd[s], center, hw)
            if st is not None:
                dens_by_hw.append(st[0])
        if not dens_by_hw:
            continue
        st_top = band_stats(psd[s], center, hw_top)   # flatness at the widest usable halfwidth
        if st_top is None:
            continue
        worst_db = st_top[1]
        records[s].append(dict(center=center,
                               density=float(np.median(dens_by_hw)),
                               worst_db=worst_db,
                               bw_ratio=max(dens_by_hw) / min(dens_by_hw),   # ~1 => bandwidth-invariant
                               clean=(worst_db <= PEAK_DB)))


# ---- report ------------------------------------------------------------------
# consensus floor per sensor = median density over bands that are BOTH flat and
# bandwidth-invariant (i.e. genuinely empty), scaled to the demod band.
summary = {}
for s in SENSORS:
    good = [r for r in records[s] if r['clean'] and r['bw_ratio'] < 2.0]
    dens = np.array(sorted(r['density'] for r in good))
    if len(dens):
        floor = float(np.median(dens))
        q1, q3 = np.percentile(dens, [25, 75])
        summary[s] = dict(floor=floor, var=floor * DEMOD_NOISE_BW,
                          rms_uv=np.sqrt(floor * DEMOD_NOISE_BW) * 1000,
                          n_good=len(good), n_rej=len(records[s]) - len(good),
                          spread=(q3 / q1 if q1 > 0 else np.inf),
                          examples=[r['center'] for r in sorted(good, key=lambda r: r['worst_db'])[:3]])
    else:
        summary[s] = None

quiet = min((d['floor'] for d in summary.values() if d), default=None)

# ---- carrier level & per-channel SNR ----------------------------------------
# Signal = the carrier tone at CARRIER_FREQ (integrate its PSD peak -> a^2/2).
# Noise  = the channel's own floor scaled to DEMOD_NOISE_BW, i.e. Var(n).
# SNR is reported in that demod band (the working band the cross-sensor cell
# uses) as carrier power / noise power. This is the raw per-channel carrier SNR,
# not the inter-diode PUF comparison SNR (that lives in the cross-sensor cell).
def carrier_power(pxx, floor_density):        # integrated tone power at the carrier, mV^2
    m = (freqs >= CARRIER_FREQ - CARRIER_HW) & (freqs <= CARRIER_FREQ + CARRIER_HW)
    if m.sum() < 3:
        return None
    return float(np.trapezoid(np.clip(pxx[m] - floor_density, 0.0, None), freqs[m]))

for s in SENSORS:
    d = summary[s]
    if d is None:
        continue
    p_car = carrier_power(psd[s], d['floor'])          # mV^2, tone power = a_i^2 / 2
    if p_car is None or p_car <= 0:
        d['a_mv'], d['snr_db'] = None, None
        continue
    d['a_mv']   = float(np.sqrt(2.0 * p_car))          # carrier amplitude a_i, mV
    d['snr_db'] = float(10.0 * np.log10(p_car / d['var']))   # 10log10( (a^2/2) / Var(n) ), demod band

# ---- readable summary (printed, and saved as summary.txt) --------------------
result  = "=" * 68 + "\n"
result += " NOISE FLOOR  —  the noise each photodiode adds on its own\n"
result += "=" * 68 + "\n"
result += "Extra electronic/pickup noise per channel, measured only where no\n"
result += "carrier or interference lives, then scaled to the analysis band.\n"
result += f"Record : {n_samp} samples ({n_samp/Fs:.0f} s)   Fs = {Fs:.0f} Hz   carrier = {CARRIER_FREQ:.1f} Hz\n"
result += f"Scaled : to the {DEMOD_NOISE_BW:.0f} Hz demod window (same band the cross-sensor cell uses)\n\n"

result += (f"  {'sensor':<9}{'carrier mV':>12}{'noise uV rms':>13}{'SNR dB':>9}"
           f"{'Var(n) mV^2':>15}{'vs quietest':>13}{'consistency':>13}\n")
result += "  " + "-" * 84 + "\n"
for s in SENSORS:
    d = summary[s]
    if d is None:
        result += f"  {s:<9}{'— no clean band found —':>54}\n"
        continue
    tag     = "consistent" if d['spread'] < 2 else "varies*"
    a_str   = f"{d['a_mv']:.3f}"   if d.get('a_mv')   is not None else "—"
    snr_str = f"{d['snr_db']:.1f}" if d.get('snr_db') is not None else "—"
    result += (f"  {s:<9}{a_str:>12}{d['rms_uv']:>13.1f}{snr_str:>9}"
               f"{d['var']:>15.2e}{d['floor']/quiet:>12.1f}x{tag:>13}\n")

result += "\nHow to read this\n"
result += "  carrier mV    carrier amplitude a_i, from integrating the PSD peak at the\n"
result += "                carrier (tone power = a_i^2/2). This is the wanted signal.\n"
result += "  noise uV rms  RMS of the channel's own noise over the demod band, in µV.\n"
result += "                It is sqrt(Var(n)) — a standard deviation (zero-mean), not a\n"
result += "                peak and not a per-√Hz density.\n"
result += "  SNR dB        carrier power over noise power in the demod band,\n"
result += "                10·log10((a_i^2/2)/Var(n)) — i.e. carrier-RMS / noise-RMS in dB.\n"
result += "                Per-channel carrier SNR, not the inter-diode PUF comparison.\n"
result += "  Var(n) mV^2   the same noise as a variance — it IS the Var(n_i) the\n"
result += "                cross-sensor cell prints for this diode.\n"
result += "  vs quietest   how many times noisier than the best channel.\n"
result += "  consistency   'consistent' = same floor everywhere; 'varies*' = the floor\n"
result += "                changes with frequency (usually a 1/f rise down low), so read\n"
result += "                the number as a typical mid-band value, not an exact constant.\n"

result += "\nWhere each number came from\n"
for s in SENSORS:
    d = summary[s]
    if d is None:
        continue
    ex = ", ".join(f"{c:.0f}" for c in d['examples'])
    result += (f"  {s:<9} averaged over {d['n_good']} clean bands "
               f"(e.g. {ex} Hz);  {d['n_rej']} rejected as too near a spur\n")
result += "  full per-band table saved alongside as summary_verbose.txt\n"

# ---- verbose per-band table (saved to its own file, never printed) -----------
verbose  = "FULL PER-BAND NOISE-FLOOR TABLE\n"
verbose += f"Fs = {Fs:.2f} Hz  carrier = {CARRIER_FREQ:.3f} Hz  Nyquist = {nyq:.1f} Hz\n"
verbose += f"record = {n_samp} samples ({n_samp/Fs:.1f} s)  Welch NFFT = {nper}  bin = {bin_hz:.3f} Hz\n"
verbose += f"candidate carrier-free bands = {len(cand)}  (guard ±{GUARD_HZ} Hz from every comb line)\n"
verbose += f"'clean' = no bin > band-median + {PEAK_DB} dB;  bw_ratio = density spread across halfwidths {BAND_HALFWIDTHS}\n"
for s in SENSORS:
    verbose += f"\n{'=' * 72}\n{s}\n{'=' * 72}\n"
    verbose += f"{'center Hz':>10}{'density mV^2/Hz':>18}{'worst dB':>10}{'bw_ratio':>10}{'clean':>8}\n"
    for r in sorted(records[s], key=lambda r: r['density']):
        verbose += (f"{r['center']:>10.1f}{r['density']:>18.3e}{r['worst_db']:>10.1f}"
                    f"{r['bw_ratio']:>10.2f}{('yes' if r['clean'] else 'NO*'):>8}\n")
    d = summary[s]
    if d:
        verbose += (f"  consensus floor = {d['floor']:.3e} mV^2/Hz over {d['n_good']} clean bands  ->  "
                    f"Var(n) = {d['var']:.4f} mV^2  (RMS {d['rms_uv']:.1f} uV)\n")
        if d.get('a_mv') is not None:
            verbose += (f"  carrier a = {d['a_mv']:.3f} mV  ->  SNR = {d['snr_db']:.1f} dB "
                        f"in the {DEMOD_NOISE_BW:.0f} Hz demod band\n")

print(result)
save_dir = save_run_data(FOLDER, data=[data], filenames=["data.csv"], sample_count=NUM_SAMPLES_TO_SAVE, summary=result)
save_summary(save_dir, verbose, "summary_verbose.txt")

if PLOT_SUMMARY:
    plot_noise_floor_summary(summary, records, DEMOD_NOISE_BW, n_samp, save_dir=save_dir)

