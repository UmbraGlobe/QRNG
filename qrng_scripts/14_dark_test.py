r"""
## Dark Test: Carrier-Band Electrical Pickup

**Purpose:** Measures deterministic electrical coupling at the carrier when no optical signal should be present. This catches PWM, SPI, or output-enable pickup that ordinary magnitude demodulation can hide.

**Physical preparation:** Keep the carrier-generating firmware running, turn all LEDs off, block ambient light, and leave the sensors in the normal electrical configuration.

**Usage:** Set `NUM_SAMPLES_TO_SAVE` and keep `CARRIER_FREQ` pinned to the known firmware carrier (currently 930.5 Hz). **Do not recalibrate on the dark record.** Set `DEMOD_BW` to the total carrier analysis width. Set `LIT_A_MV` to a representative lit carrier amplitude, or `None` to omit that comparison. Then run with live data or the first saved dark file.

`pickup uV` is coherent phase-locked amplitude; random noise averages out. Similar phase across channels indicates a shared electrical source. Compare pickup with `in-band uV` and the lit carrier because pickup above the independent floor can bias cross-sensor covariance. Live mode saves raw data, `summary.txt`, and a figure under `data/dark_test_data/...`.

Source: notebook cell 31 (code copied unchanged unless noted below).
"""

from qrng_common import *
from qrng_viz import *

### CONFIG

CARRIER_FREQ       = 930.5
NUM_SAMPLES_TO_SAVE = 500000  
DEMOD_BW            = 4.0   
FOLDER             = "dark_test_data"
LIT_A_MV           = 1052.6   # mean a_i from the last LIT run, for the "% of carrier" tile (None to skip)
PLOT_SUMMARY       = True     # render the summary figure (saved next to summary.txt)

###

halfbw = DEMOD_BW / 2.0
data   = read_from_file_instant(LOAD_FILE_PATHS[0]) if LOAD_FROM_FILE else read_data_realtime(NUM_SAMPLES_TO_SAVE)
data   = data[:NUM_SAMPLES_TO_SAVE]
n_samp = len(data)

n      = np.arange(n_samp)
lo     = np.exp(-1j * 2 * np.pi * CARRIER_FREQ * n / Fs)
sos_lp = signal.butter(4, halfbw, 'low', fs=Fs, output='sos')

pick = {}
for s in SENSORS:
    v  = adc_to_mv(data[s].astype(float))
    v  = v - v.mean()                             # kill the DC pedestal so it can't leak into c_i
    z  = v * lo
    c  = 2.0 * np.mean(z)                          # coherent carrier-band phasor c_i (mV)
    zb = signal.sosfiltfilt(sos_lp, z)
    pick[s] = dict(c=c, mag=abs(c), phase=np.degrees(np.angle(c)),
                   noise=2.0 * np.std(zb - zb.mean()))

ref = SENSORS[0]

result  = "=" * 70 + "\n"
result += " DARK TEST  —  carrier-band electrical pickup (LEDs OFF)\n"
result += "=" * 70 + "\n"
result += "Coherent phasor at the carrier with no light: what couples in from the carrier\n"
result += "PWM / SPI / OE switching. Complex baseband (a quadrature spur is invisible to\n"
result += "magnitude demod).\n"
result += f"Record : {n_samp} samples ({n_samp/Fs:.0f} s)   Fs = {Fs:.0f} Hz   carrier = {CARRIER_FREQ:.1f} Hz\n\n"

hdr = f"  {'sensor':<11}{'pickup uV':>11}{'phase deg':>11}{'in-band uV':>12}"
result += hdr + "\n  " + "-" * (len(hdr) - 2) + "\n"
for s in SENSORS:
    d = pick[s]
    line = f"  {s:<11}{d['mag']*1000:>11.2f}{d['phase']:>11.1f}{d['noise']*1000:>12.2f}"
    result += line + "\n"

result += f"\nCommon-mode check (phasor vs {ref}):\n"
for s in SENSORS:
    dphi = np.degrees(np.angle(pick[s]['c'] / pick[ref]['c'])) if pick[ref]['mag'] > 0 else float('nan')
    result += f"  {s:<11} |c| = {pick[s]['mag']*1000:8.2f} uV   dphase = {dphi:+7.1f} deg\n"
result += "  aligned dphase across channels => one shared source, phase-locked to the carrier.\n"

result += "\nVerdict\n"
worst = max(pick.values(), key=lambda d: d['mag'])
result += f"  largest pickup = {worst['mag']*1000:.2f} uV.\n"
result += "  Compare pickup uV to the ~80-210 uV noise floor: below it = buried; above = a coherent\n"
result += "  line that WILL bias Var(m)/cov in the cross-sensor cell (subtract it, per analysis.md).\n"

result += "\nHow to read this\n"
result += "  pickup uV    |c_i| = 2*mean(v*e^-jwt): deterministic carrier-band pickup amplitude.\n"
result += "               Random noise averages out; only phase-locked coupling survives.\n"
result += "  phase deg    arg(c_i) vs the demod LO. Compare to the LIT carrier phase (variance cell\n"
result += "               de-rotation ref) to tell an in-phase pickup from a quadrature one.\n"
result += f"  in-band uV   incoherent fluctuation in the +/-{halfbw:.0f} Hz demod band (2*std of baseband).\n"
result += "  % of carrier |c_i| relative to the lit-run a_i -- how much it perturbs that channel.\n"

print(result)
save_dir = save_run_data(FOLDER, data=[data], filenames=["data.csv"], sample_count=NUM_SAMPLES_TO_SAVE, summary=result)   # CHANGED: keep the run folder

if PLOT_SUMMARY:
    plot_dark_test_summary(pick, n_samp, halfbw, lit_a=LIT_A_MV, save_dir=save_dir)
