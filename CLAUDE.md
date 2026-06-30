# QRNG — Photodiode Noise Characterization

Experimental project working toward a **quantum random number generator (QRNG)** based on
photodiode noise. The immediate goal is *characterization*: measure each photodiode's
fluctuations and separate the **independent per-diode noise** (shot/thermal — the entropy
source we want) from the **common-mode classical LED intensity fluctuation** (correlated,
predictable — the thing we must reject).

All work lives in [main.ipynb](main.ipynb).

## Hardware & data pipeline

- **MCU**: Arduino streaming over serial — `PORT = "COM8"`, `BAUD = 230400`.
- **Sensors**: 4 photodiodes — `['BPW34_1', 'BPW34_2', 'Inano_1', 'Inano_2']`. BPW34 are
  silicon PIN diodes; "Inano" are a different diode type. Both pairs view the **same LED**.
- **LED**: driven at a **500 Hz carrier** (`CARRIER_FREQ`). It is *amplitude-modulated*, not
  DC — the 500 Hz tone is the carrier; intensity fluctuations ride on it as sidebands. Running
  at 500 Hz (vs. DC) deliberately escapes 1/f flicker noise and drift near baseband.
- **ADC**: 14-bit, `ADC_MAX = 16383`, `V_REF = 5`. Convert with `adc_to_mv()`.
- **Sample rate**: ~1500 Hz (`Fs = 1500` hardcoded; measured ≈ 1499.25 Hz from timestamps).
- **Wire format**: packets of `0xAA 0xBB` sync bytes + `SAMPLES_PER_PACKET = 50` records of
  `payload_dtype` = `(timestamp u4, BPW34_1 u2, BPW34_2 u2, Inano_1 u2, Inano_2 u2)`.
- **Storage**: CSV under `data/<MM_DD_YYYY>/run_<HH_MM_SS>.csv`; multi-window captures under
  `VarianceData/<date>/run_<time>/sampleN.csv`.

### Read paths (set `LOAD_FROM_FILE`)
- `try_read_realtime(ser)` — read one 50-sample packet from serial.
- `read_from_file_realtime(path)` — replay a CSV at wall-clock pace (for the live plots).
- `read_from_file_instant(path)` — load a whole CSV at once (for batch analysis).
- `save_data(n, path, ser)` — capture `n` samples to CSV.

## Notebook cell map

1. **Setup** — imports, constants (UPPER_CASE), `payload_dtype`, helper functions. Run first.
2. **Saves realtime data to file** — capture to CSV.
3. **Realtime plot** — scrolling pyqtgraph view of the 4 raw channels (mV).
4. **Spectrogram Plotting** — per-channel spectrograms; confirm the carrier sits at 500 Hz.
5. **Realtime PSD** — live Welch PSD per channel, optional bandpass.
6. **Find the Variance Between Sensors** — the core analysis (model + covariance method).
7. **Variance of A Single Sensor over Sequential Time Periods** — stationarity / repeatability
   check on one sensor across time windows.

## The model & math (cell 6)

Measured load voltage for diode *i*:

```
V_i(t) = A_i(t) · cos(2π f_c t + θ)
A_i(t) = a_i + a_i·m(t) + n_i(t)
```

- `a_i` — **constant** mean carrier amplitude (mV). (It is a mean, so not time-varying.)
- `m(t)` — LED **relative** intensity fluctuation, dimensionless, `E[m] = 0`, **common to all
  diodes** with coupling proportional to `a_i` (each sees the *same fractional* fluctuation).
- `n_i(t)` — independent in-band noise of each diode (shot, thermal), mV, uncorrelated across
  diodes and with `m`.

**Recover `A_i(t)`** by I/Q (lock-in) demodulation of the 500 Hz carrier, then magnitude:
`A = 2·√(LP[V·cos]² + LP[V·sin]²)`. The ×2 makes the envelope equal the carrier amplitude.

**Separate common vs. independent noise** via cross-covariance (independent noise cancels):

```
Cov(A_i, A_j) = a_i · a_j · Var(m)          ⟹  Var(m) = Cov(A_i,A_j) / (a_i·a_j)
Var(A_i)      = a_i² · Var(m) + Var(n_i)     ⟹  Var(n_i) = Var(A_i) − a_i²·Var(m)
```

`Var(m)` is estimated as the **mean of the normalized cross-covariances over correlated pairs
only** (`|corr| > CORR_GATE`, default 0.5) so a decoupled channel can't bias it.

Three matrices are printed: `[1]` covariance `cov_A` (mV²; diag = total per-channel variance),
`[2]` correlation `corr_A`, `[3]` normalized covariance `cov/(a_i a_j)` (off-diagonals would be
a constant `= Var(m)` if a single fractional `m` explained everything).

## Key empirical findings

- **The BPW34 pair is almost entirely common-mode.** `corr(BPW34_1, BPW34_2) ≈ 0.97`; their
  independent noise is only ~1–9% of total variance. ~90%+ of what they measure is the shared
  LED/driver fluctuation.
- **Inano_2 is effectively decoupled** from the LED common mode (near-zero / slightly negative
  correlation with the others, tiny carrier amplitude ~0.5 mV). Its `% noise` is flagged
  unreliable because the `a_i²·Var(m)` subtraction assumes coupling = amplitude, which fails.
- **`Var(m)` RMS ≈ 0.25 (25%).** This is *far* too large to be intrinsic LED RIN or shot noise
  (real LED RIN ≪ 1%). It is a **classical** fluctuation — slow brightness drift (warm-up /
  thermal) and/or an unregulated/PWM driver and supply ripple. With a ±2 Hz bandpass, slow
  drift passes straight into the envelope and dominates.
- **Implication for QRNG**: the classical common mode beats the private noise by ~5× in
  amplitude (~25× in power), so a single raw channel is unusable as entropy. The intended fix is
  **common-mode rejection** — e.g. `BPW34_1/a_1 − BPW34_2/a_2` cancels `m(t)` and leaves the
  combined private noise `Var(n_1)/a_1² + Var(n_2)/a_2²`, with no large-number subtraction.
- **Stationarity (cell 7)**: same sensor over 6 sequential windows gives `Var(A)` spread
  (CV) ≈ 5% — the total fluctuation is reasonably stationary window-to-window. Note this cell
  measures **one** detector, so it *cannot* separate `m` from `n` (no second channel to cancel
  noise); it's a repeatability check only.

## Gotchas / lessons learned (don't re-introduce these)

- **`filtered_all[s]` is a list of per-packet chunks**, not one array. Concatenate first, *then*
  slice samples: `np.concatenate(filtered_all[s])[200:]` — slicing the list (`[200:]`) drops
  whole packets and empties short runs.
- **`fir_lowpass_filter` returns `(filtered, taps)`** — unpack it (`I_f, _ = ...`). Using the
  tuple directly breaks `I_f**2`. The helper wraps `lfilter` in `np.asarray` so the return type
  is unambiguously `ndarray` (fixes a Pylance `reportOperatorIssue`).
- **Filter transients**: both the Butterworth bandpass (`sosfilt` startup) and the 101-tap FIR
  ring at the start. Drop ~200 samples post-bandpass and ~150 post-FIR before computing stats.
- **Magnitude demod is biased** (Rician): `2√(I²+Q²)` rectifies noise, so single-channel
  envelope variance has a noise floor and `a_i` is slightly inflated. This largely *cancels* in
  the cross-covariance (noise is independent), so trust the cross terms, not a lone channel.
- **Negative `Var(n)`** (seen on BPW34_2) is unphysical and means `a_i²·Var(m) > Var(A_i)` — a
  catastrophic-cancellation artifact for a channel that is ~100% common-mode (subtracting two
  nearly-equal large numbers). It signals "independent noise below the method's resolution,"
  i.e. consistent with ~0. Use the matched-pair `Var(m)` or the difference method instead, and/or
  longer averaging.
- **Very narrow bandpass is risky**: an 8th-order Butterworth at 500±0.5 Hz has poles almost on
  the unit circle (numerically touchy, very long ringing). ±2 Hz is the current compromise.
- **`Var(m)` averaging must gate on correlation** — including a decoupled channel drags the mean
  down. `CORR_GATE` controls this.
- **`Fs = 1500` is hardcoded** in the analysis cells though the true rate (~1499.25 Hz) is
  measured from timestamps. The magnitude demod is robust to small `Fs` error, but the bandpass
  centering and `t` axis are not — keep this in mind if the carrier looks off-center.

## Conventions

- Constants in `UPPER_CASE` (defined in the Setup cell); helper functions also live there.
- Envelope/analysis dicts are keyed by sensor name (or window index in cell 7).
- `var(m)` is **dimensionless** (mV²/mV²); its RMS is the fractional intensity fluctuation.
- `Var(A)`, `Var(n)`, `cov_A` are in mV².

## Open TODOs (from in-notebook comments)

1. Confirm `V_REF`; check whether samples-per-carrier-cycle matters; confirm the LED drive is a
   clean sine.
2. Support read-from-file for every analysis block.
3. Cross-check by sampling multiple blocks from the **same** PD at different times (cell 7) vs.
   the across-PD method (cell 6).
4. Sweep bandpass bands.
5. Determine whether the ~25% RMS common mode is slow drift (fixable by AC-coupling / high-pass /
   regulated current source) or genuine driver noise — high-pass the envelope and see if RMS
   collapses.
6. Implement the common-mode-rejected difference signal as the entropy candidate and test whether
   its noise scales with optical power (shot-noise-like) vs. fixed (electronic).
