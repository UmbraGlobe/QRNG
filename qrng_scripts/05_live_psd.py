r"""
## Live Power Spectrum

**Purpose:** Continuously estimates the power spectral density (PSD) of every channel to locate the carrier, harmonics, packet-related lines, and broadband noise.

**Usage:** Set `SAMPLE_MEMORY` to the rolling history size. Keep `FILTER = False` to see the full spectrum from 0 Hz to Nyquist; use `True` only to inspect the carrier band selected by `BANDWIDTH`. Run the cell and close the window to stop it.

The vertical axis is logarithmic in `mV²/Hz`: narrow tall peaks are tones and the lower background is the noise floor. Saved-file mode replays the first file and stops at end of file. This on-screen diagnostic does not save results.

Source: notebook cell 13 (code copied unchanged unless noted below).
"""

from qrng_common import *

### CONFIG

FILTER = False # set to True to enable bandpass filtering
SAMPLE_MEMORY = 2000
BANDWIDTH = 4  # Hz, width of the bandpass filter around the carrier frequency

###

buffers = {s: np.zeros(SAMPLE_MEMORY) for s in SENSORS}
if FILTER:
    sos_bp = signal.butter(4, [CARRIER_FREQ - BANDWIDTH/2, CARRIER_FREQ + BANDWIDTH/2], btype='bp', fs=Fs, output='sos')
    sos_zi = {s: signal.sosfilt_zi(sos_bp) for s in SENSORS}

app = pg.mkQApp("QRNG PSD realtime")
win = pg.GraphicsLayoutWidget(show=True, title="QRNG PSD realtime")
curves = {}
plots = []
for i, s in enumerate(SENSORS):
    p = win.ci.addPlot(row=0, col=i, title=s)
    p.setLabel('left', 'PSD (mV^2/Hz)')
    p.setLabel('bottom', 'Frequency (Hz)')
    p.setLogMode(x=False, y=True)              # log Y: carrier + floor both visible, DC no longer blows up the scale
    p.getViewBox().setXRange(0, Fs / 2, padding=0)   # 0 -> Nyquist, tracks whatever Fs is

    plots.append(p)
    curves[s] = p.plot()

if LOAD_FROM_FILE:
    reader = read_from_file_realtime(LOAD_FILE_PATHS[0])

def update():
    if LOAD_FROM_FILE:
        data = reader()
    else:
        data = try_read_realtime(ser)
    if data is None:
        if LOAD_FROM_FILE and reader.exhausted:
            timer.stop()
        return
    if data is not None:
        for s in SENSORS:
            chunk_mv = adc_to_mv(data[s])
            buffers[s] = np.roll(buffers[s], -SAMPLES_PER_PACKET)
            if FILTER:
                chunk, sos_zi[s] = signal.sosfilt(sos_bp, chunk_mv, zi=sos_zi[s])
                buffers[s][-SAMPLES_PER_PACKET:] = chunk
            else:
                buffers[s][-SAMPLES_PER_PACKET:] = chunk_mv
            f, pxx = welch(buffers[s], fs=Fs, nperseg=512)
            curves[s].setData(f, pxx)

timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start()
try:
    with serial_port() as ser:          # Fix 1: port always closed, even on window-close/interrupt
        if not LOAD_FROM_FILE:
            drain_to_latest(ser)
        pg.exec()
finally:
    timer.stop()
