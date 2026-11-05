r"""
## Live Envelope Ripple Monitor

**Purpose:** Demodulates every channel and displays slow amplitude ripple on the carrier envelope instead of the raw carrier waveform.

**Usage:** Normally keep `DEMOD_FREQ = CARRIER_FREQ`. Set `ENV_LP_HZ` high enough to pass the ripple frequencies of interest. Set `MEASURE_AT_HZ` to a known frequency for a calibrated readout, or leave it `None` to report the strongest non-DC ripple tone. `SAMPLE_MEMORY` must be a multiple of 50. Run the cell and close the window to stop it.

Each title reports mean envelope amplitude, total RMS ripple as a percentage, the strongest ripple tone, and the optional requested tone. Wait for the buffer and filters to settle before reading values. Saved-file mode replays the first file. This diagnostic does not save results.

Source: notebook cell 29 (code copied unchanged unless noted below).
"""

from qrng_common import *

### CONFIG

DEMOD_FREQ    = CARRIER_FREQ   # frequency to lock in on and measure ripple at (Hz). Try 625 vs CARRIER_FREQ.
ENV_LP_HZ     = 50.0           # I/Q low-pass: passes envelope ripple below this, kills the 2*f0 sum image.
MEASURE_AT_HZ = None           # report ripple amplitude at this exact freq too (None = auto-pick the peak).
SAMPLE_MEMORY = 10000          # envelope samples shown (multiple of SAMPLES_PER_PACKET)
STATS_EVERY   = 10             # refresh the ripple readout every N packets (throttle so the label doesn't flicker)

###

assert SAMPLE_MEMORY % SAMPLES_PER_PACKET == 0

# Lock-in envelope path: A = 2*sqrt(LP[v*cos]^2 + LP[v*sin]^2), matching the analysis cell.
# The LP removes the 2*DEMOD_FREQ image; anything left below ENV_LP_HZ is genuine envelope ripple.
sos_lp = signal.butter(4, ENV_LP_HZ, btype='low', fs=Fs, output='sos')
zi_I = {s: np.zeros_like(signal.sosfilt_zi(sos_lp)) for s in SENSORS}
zi_Q = {s: np.zeros_like(signal.sosfilt_zi(sos_lp)) for s in SENSORS}

env = {s: np.full(SAMPLE_MEMORY, np.nan) for s in SENSORS}   # demodulated envelope A(t), mV
write_idx = 0
demod_n   = 0        # running sample index -> keeps the LO phase continuous across packet boundaries
pkt_count = 0
label_style = {'color': '#FFF', 'font-size': '14pt'}

def ripple_stats(a, target_hz):
    """RMS ripple (% of mean), the dominant ripple tone, and (optionally) the tone at target_hz,
    over a settled envelope block. Hann + coherent-gain correction for calibrated amplitudes."""
    a = a[np.isfinite(a)]
    if a.size < 512 or np.mean(a) == 0:
        return None
    mean = np.mean(a)
    ac   = a - mean
    pct  = 100.0 * np.std(ac) / abs(mean)                 # total RMS ripple, % of the mean envelope
    win  = np.hanning(a.size)
    Xw   = np.abs(np.fft.rfft(ac * win))
    fb   = np.fft.rfftfreq(a.size, d=1.0 / Fs)
    cg   = np.sum(win)                                    # coherent gain: amp = 2|X_k| / cg
    kpk  = 1 + int(np.argmax(Xw[1:]))                     # skip the DC bin
    peak_f, peak_amp = fb[kpk], 2.0 * Xw[kpk] / cg
    tgt = None
    if target_hz is not None:
        kt  = int(np.argmin(np.abs(fb - target_hz)))
        tgt = (fb[kt], 2.0 * Xw[kt] / cg)
    return mean, pct, peak_f, peak_amp, tgt

app = pg.mkQApp("ripple")
win_gl = pg.GraphicsLayoutWidget(show=True, title="ripple")
pen = pg.mkPen('r', width=2)
curves, plots = {}, {}
for i, s in enumerate(SENSORS):
    p = win_gl.ci.addPlot(row=i, col=0)
    p.setLabel('left', 'A (mV)', **label_style)
    p.setLabel('top', f'{s}', **label_style)
    plots[s] = p
    curves[s] = p.plot(pen=pen, connect='finite')  # 'finite' leaves a gap at the nan cursor

if LOAD_FROM_FILE:
    reader = read_from_file_realtime(LOAD_FILE_PATHS[0])

def update():
    # Sweep display of the demodulated envelope: each packet overwrites 50 envelope samples at the cursor.
    global write_idx, demod_n, pkt_count
    got_packet = False
    while True:
        if LOAD_FROM_FILE:
            data = reader()
            if data is None:
                if reader.exhausted:
                    timer.stop()
                break
        else:
            if ser is not None and ser.in_waiting < BYTES_PER_SAMPLE + 2:
                break
            data = try_read_realtime(ser)
            if data is None:
                continue
        # continuous LO phase across packets -- demod_n never resets mid-run (a phase jump corrupts I/Q)
        idx  = demod_n + np.arange(SAMPLES_PER_PACKET)
        lo_c = np.cos(2 * np.pi * DEMOD_FREQ * idx / Fs)
        lo_s = np.sin(2 * np.pi * DEMOD_FREQ * idx / Fs)
        for s in SENSORS:
            v = adc_to_mv(data[s])
            I, zi_I[s] = signal.sosfilt(sos_lp, v * lo_c, zi=zi_I[s])
            Q, zi_Q[s] = signal.sosfilt(sos_lp, v * lo_s, zi=zi_Q[s])
            env[s][write_idx:write_idx + SAMPLES_PER_PACKET] = 2.0 * np.sqrt(I**2 + Q**2)
        demod_n  += SAMPLES_PER_PACKET
        write_idx = (write_idx + SAMPLES_PER_PACKET) % SAMPLE_MEMORY
        got_packet = True
        if LOAD_FROM_FILE:
            break
    if got_packet:
        pkt_count += 1
        show_stats = (pkt_count % STATS_EVERY == 0)
        for s in SENSORS:
            env[s][write_idx:write_idx + SAMPLES_PER_PACKET] = np.nan  # blank slot marks the cursor
            curves[s].setData(env[s])
            if show_stats:
                st = ripple_stats(env[s], MEASURE_AT_HZ)
                if st is not None:
                    mean, pct, pf, pa, tgt = st
                    txt = f'{s}   A={mean:.3f} mV   ripple={pct:.2f}%   peak {pf:.1f} Hz ({pa:.4f} mV)'
                    if tgt is not None:
                        txt += f'   @{tgt[0]:.1f} Hz: {tgt[1]:.4f} mV'
                    plots[s].setLabel('top', txt, **label_style)

timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(10)
try:
    with serial_port() as ser:          # port always closed, even on window-close/interrupt
        if not LOAD_FROM_FILE:
            drain_to_latest(ser)
        pg.exec()
finally:
    timer.stop()  # don't leave a stale timer attached if the cell is re-run
