r"""
## Carrier and Aliasing Safety Check

**Purpose:** Predicts where PWM-carrier harmonics fold after sampling and checks whether a significant harmonic lands inside the demodulation band.

**Usage:** Run after calibration and whenever sample rate, carrier, PWM duty, or analysis bandwidth changes. `F0` and `FS` normally inherit `CARRIER_FREQ` and `Fs`. Set `PWM_DUTY` to the firmware duty cycle and `GUARD_HZ` to the demodulation half-width. A 50% duty cycle suppresses even harmonics.

`[OK]` means the fundamental is below Nyquist and no modeled significant harmonic aliases into the guarded carrier band. `[FAIL]` means the timing plan can contaminate demodulation; correct it before trusting later analysis. This is an analytical prediction and does not capture or save data.

Source: notebook cell 25 (code copied unchanged unless noted below).
"""

from qrng_common import *

### CONFIG

F0        = CARRIER_FREQ   # carrier fundamental (Hz), estimated in Setup
FS        = Fs             # sample rate (Hz), estimated in Setup
PWM_DUTY  = 0.50           # 0.50 => even harmonics vanish (phase2 ruling: keep 50%)
GUARD_HZ  = 2.0            # demod half-band around the carrier (±2 Hz working band)
DBC_FLOOR = -40.0          # ignore harmonics weaker than this vs. the fundamental
MAX_ORDER = 4000           # sweep depth (comb is endless; DBC_FLOOR does the culling)

###

nyq = FS / 2.0

def fold(f, fs):
    """Alias a real tone at frequency f into [0, fs/2]."""
    f = np.abs(f) % fs
    return fs - f if f > fs / 2 else f

def pwm_dbc(k, duty):
    """k-th harmonic of a duty-D pulse train, in dBc vs. the fundamental.
    np.sinc is the normalized sinc: |a_k| ∝ |sinc(k*duty)|."""
    ak, a1 = abs(np.sinc(k * duty)), abs(np.sinc(duty))
    return -np.inf if ak == 0 else 20 * np.log10(ak / a1)

# --- sweep the comb (harmonics k >= 2) ----------------------------------------
rows = []
for k in range(2, MAX_ORDER + 1):
    dbc = pwm_dbc(k, PWM_DUTY)
    if dbc < DBC_FLOOR:          # too weak to matter
        continue
    alias = fold(k * F0, FS)
    rows.append((k, k * F0, alias, dbc, abs(alias - F0)))

# --- verdicts -----------------------------------------------------------------
print(f"Fs = {FS:.3f} Hz   Nyquist = {nyq:.3f} Hz   f0 = {F0:.3f} Hz   duty = {PWM_DUTY:.0%}")
print(f"guard = ±{GUARD_HZ} Hz   dBc floor = {DBC_FLOOR}   orders swept = 2..{MAX_ORDER}\n")

if F0 < nyq:
    print(f"[OK]   fundamental below Nyquist  ({F0:.3f} < {nyq:.3f} Hz)")
else:
    print(f"[FAIL] fundamental {F0:.3f} Hz is at/above Nyquist {nyq:.3f} Hz — undersampled!")

collisions = [r for r in rows if r[4] <= GUARD_HZ]
if collisions:
    print(f"\n[FAIL] {len(collisions)} harmonic(s) alias into the ±{GUARD_HZ} Hz carrier band:")
    for k, fk, al, dbc, d in collisions:
        print(f"   k={k:<5d} {fk:11.1f} Hz -> {al:8.3f} Hz  ({dbc:+.1f} dBc, {d:.3f} Hz from f0)")
else:
    k, fk, al, dbc, d = min(rows, key=lambda r: r[4])
    print(f"\n[OK]   no harmonic within ±{GUARD_HZ} Hz of the carrier.")
    print(f"       nearest approach: k={k} ({dbc:+.1f} dBc) folds to {al:.3f} Hz, {d:.3f} Hz away.")

# --- spur table ---------------------------------------------------------------
print("\nspur table (significant harmonics):")
print(f"   {'k':>5}  {'harmonic Hz':>12}  {'alias Hz':>10}  {'dBc':>7}")
for k, fk, al, dbc, d in rows:
    print(f"   {k:5d}  {fk:12.1f}  {al:10.3f}  {dbc:+7.1f}")
