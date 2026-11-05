r"""
## Live Pairwise Bit Display

**Purpose:** Displays the latest pairwise sensor-comparison bits as a black-and-white grid. Black is 1 (the first sensor is larger), white is 0, and grey cells are unused.

**Usage:** Set `AVERAGE_SAMPLES` for coherent averaging per read. `UPDATE_MS` controls refresh speed; `CELL_PX`, `QUIET_ZONE`, and `SHOW_BITS_TEXT` only affect appearance. Run the cell and close the window to stop it and release the serial port. Saved-file mode replays the first file.

The title shows total reads, the fraction belonging to the most common bitstring, and estimated LO frequency error. The grid is QR-like visualization only, not a scannable standards-compliant QR code. Startup behavior is included in its stability count and nothing is saved; use Fingerprint Averaging and Stability for a report.

Source: notebook cell 23 (code copied unchanged unless noted below).
"""

from qrng_common import *

### CONFIG

AVERAGE_SAMPLES  = 1       # samples coherently averaged per fingerprint read (1 = a read per sample)
DEMOD_LP_HZ      = 20       # Hz, low-pass cutoff for the demodulated envelope
CELL_PX          = 60       # on-screen size of one code cell -> the window scales with the bitstring
QUIET_ZONE       = 1        # cells of white border around the code, as a real QR has
SHOW_BITS_TEXT   = True     # digits + a_hat under the code (auto-hidden once the bitstring is long)
UPDATE_MS        = 100

###

# Adaptive geometry: L sensors -> L*(L-1)/2 comparison bits, laid out in the squarest grid that
# holds them. Nothing below is hardcoded to 4 channels, so adding photodiodes just grows the code.
PAIRS  = [(i, j) for i in range(len(SENSORS)) for j in range(i + 1, len(SENSORS))]  # same order as the averaging cell
N_BITS = len(PAIRS)
SIDE   = int(np.ceil(np.sqrt(N_BITS)))                 # columns
ROWS   = int(np.ceil(N_BITS / SIDE))

print(f"{len(SENSORS)} sensors -> {N_BITS} bits in a {ROWS}x{SIDE} grid   "
      f"{AVERAGE_SAMPLES} samples/read ({AVERAGE_SAMPLES/Fs:.4f} s, {Fs/AVERAGE_SAMPLES:.1f} reads/s)")

sos_lps = {s: signal.butter(4, DEMOD_LP_HZ, btype='low', fs=Fs, output='sos') for s in SENSORS}
zi_I = {s: np.zeros((len(sos_lps[s]), 2)) for s in SENSORS}  # type: ignore
zi_Q = {s: np.zeros((len(sos_lps[s]), 2)) for s in SENSORS}  # type: ignore
sample_count = 0


def filter_pipeline(data):
    # Returns the COMPLEX baseband Z = I + jQ (identical to the averaging cell).
    global sample_count
    N = SAMPLES_PER_PACKET
    t = np.arange(sample_count, sample_count + N) / Fs
    lo_cos, lo_sin = np.cos(2*np.pi*CARRIER_FREQ*t), -np.sin(2*np.pi*CARRIER_FREQ*t)
    sample_count += N          # the pipeline owns the counter - never advance it again in the loop

    Z = {}
    for s in SENSORS:
        v_mv = adc_to_mv(data[s])
        I, zi_I[s] = signal.sosfilt(sos_lps[s], lo_cos * v_mv, zi=zi_I[s])
        Q, zi_Q[s] = signal.sosfilt(sos_lps[s], lo_sin * v_mv, zi=zi_Q[s])
        Z[s] = I + 1j * Q
    return Z


# --- display: the bitstring as a QR-style code, black = 1, white = 0 ------------------
WHITE, BLACK, UNUSED = 0, 1, 2                                   # UNUSED = grid cells past the last bit
LUT = np.array([[255, 255, 255], [0, 0, 0], [150, 150, 150]], dtype=np.ubyte)
IMG_SHAPE = (ROWS + 2*QUIET_ZONE, SIDE + 2*QUIET_ZONE)
label_style = {'color': '#FFF', 'font-size': '16pt'}

app = pg.mkQApp("fingerprint realtime")
win = pg.GraphicsLayoutWidget(show=True, title="fingerprint realtime")
win.resize(max(420, IMG_SHAPE[1] * CELL_PX), max(420, IMG_SHAPE[0] * CELL_PX + 140))

p = win.ci.addPlot(row=0, col=0)
p.setLabel('top', 'waiting for data...', **label_style)
p.hideAxis('left'); p.hideAxis('bottom')
p.setAspectLocked(True)
p.invertY(True)                                # row 0 at the top, read like a QR code
p.setMouseEnabled(x=False, y=False)
p.setXRange(0, IMG_SHAPE[1], padding=0); p.setYRange(0, IMG_SHAPE[0], padding=0)

img = pg.ImageItem()
img.setOpts(axisOrder='row-major')             # grid[row, col], not pyqtgraph's default [x, y]
img.setLookupTable(LUT)
img.setImage(np.full(IMG_SHAPE, WHITE, dtype=np.ubyte), levels=(0, 2))
p.addItem(img)
txt = win.addLabel("", row=1, col=0, size='13pt', color='#FFF')

pend       = {s: np.empty(0, dtype=complex) for s in SENSORS}   # baseband awaiting a full read
n_reads    = 0
seen       = {}         # bitstring tuple -> count, for the live stability readout
REF        = None       # amplitude reference channel, fixed on the first read for phase continuity
last_phase = None

if LOAD_FROM_FILE:
    reader = read_from_file_realtime(LOAD_FILE_PATHS[0])


def render(bits, a_hat, stability, df):
    grid = np.full(IMG_SHAPE, WHITE, dtype=np.ubyte)
    flat = np.full(ROWS * SIDE, UNUSED, dtype=np.ubyte)         # spare cells stay grey
    flat[:N_BITS] = bits
    grid[QUIET_ZONE:QUIET_ZONE + ROWS, QUIET_ZONE:QUIET_ZONE + SIDE] = flat.reshape(ROWS, SIDE)
    img.setImage(grid, levels=(0, 2))

    p.setLabel('top', f"read {n_reads}   stable {100*stability:5.1f}%   df = {df:+.3f} Hz", **label_style)
    if SHOW_BITS_TEXT and N_BITS <= 64:
        mv = "   ".join(f"{s.split('_')[-1]}={a_hat[s]:6.2f}" for s in SENSORS)
        txt.setText(f"{''.join(str(int(b)) for b in bits)}<br>{mv} mV")


def update():
    global n_reads, REF, last_phase

    while True:                                  # drain everything buffered, same as the realtime plot cell
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
        z = filter_pipeline(data)
        for s in SENSORS:
            pend[s] = np.concatenate([pend[s], z[s]])
        if LOAD_FROM_FILE:
            break

    k = len(pend[SENSORS[0]]) // AVERAGE_SAMPLES     # complete reads available this tick
    if k == 0:
        return

    # Coherent (complex) window mean, not mean-of-magnitudes: it avoids the Rician bias. A residual
    # LO frequency error costs a sinc factor identical on every channel, so the comparisons - and
    # therefore the bits - are unaffected by it. All k reads are done at once (AVERAGE_SAMPLES = 1
    # produces one read per sample, i.e. ~Fs of them per second).
    zbar, A = {}, np.empty((k, len(SENSORS)))
    for idx, s in enumerate(SENSORS):
        zbar[s] = pend[s][:k * AVERAGE_SAMPLES].reshape(k, AVERAGE_SAMPLES).mean(axis=1)
        A[:, idx] = 2 * np.abs(zbar[s])
        pend[s] = pend[s][k * AVERAGE_SAMPLES:]

    bits_arr = np.stack([(A[:, i] > A[:, j]).astype(np.ubyte) for i, j in PAIRS], axis=1)   # (k, N_BITS)
    uniq, counts = np.unique(bits_arr, axis=0, return_counts=True)
    for row, c in zip(uniq, counts):
        key = tuple(int(b) for b in row)
        seen[key] = seen.get(key, 0) + int(c)
    n_reads += k

    if REF is None:
        REF = SENSORS[int(A[-1].argmax())]
    phs  = np.angle(zbar[REF])                   # phase walk between reads -> LO error,
    prev = phs[-2] if k > 1 else last_phase      # unambiguous to +-Fs/(2*AVERAGE_SAMPLES)
    df   = 0.0 if prev is None else \
        ((phs[-1] - prev + np.pi) % (2*np.pi) - np.pi) / (2*np.pi*AVERAGE_SAMPLES/Fs)
    last_phase = phs[-1]

    render(bits_arr[-1], {s: A[-1, i] for i, s in enumerate(SENSORS)},
           max(seen.values()) / n_reads, df)     # redraw once per tick, not once per read


timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(UPDATE_MS)
try:
    with serial_port() as ser:
        if not LOAD_FROM_FILE:
            drain_to_latest(ser)
        pg.exec()
finally:
    timer.stop()   # don't leave a stale timer attached if the cell is re-run
