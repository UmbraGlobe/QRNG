r"""
## Live Waveform Monitor

**Purpose:** Displays the latest voltage samples from all four channels so you can check connections, signal level, clipping, and obvious interference.

**Usage:** Set `SAMPLE_MEMORY` to the displayed history; it must be a multiple of 50. Keep `FILTER = False` to inspect raw signals, or use `True` to display only the carrier band selected by `BANDWIDTH`. Run the cell and close the plot window to stop it and release the serial port.

Live mode follows the Arduino. Saved-file mode replays the first file at its recorded speed and stops at end of file. The moving blank slot is the write cursor, not missing data. This diagnostic does not save a capture.

Source: notebook cell 9 (code copied unchanged unless noted below).
"""

from qrng_common import *


### CONFIG

SAMPLE_MEMORY = 10000
FILTER = True # set to True to enable bandpass filtering
BANDWIDTH = 4  # Hz, width of the bandpass filter around the carrier frequency

###

assert SAMPLE_MEMORY % SAMPLES_PER_PACKET == 0
buffers = {s: np.full(SAMPLE_MEMORY, np.nan) for s in SENSORS} 
write_idx = 0
label_style = {'color': '#FFF', 'font-size': '16pt'}

app = pg.mkQApp("realtime")
win = pg.GraphicsLayoutWidget(show=True, title="realtime")
pen = pg.mkPen('r', width=5)
curves = {}
for i, s in enumerate(SENSORS):
    p = win.ci.addPlot(row=i, col=0)
    p.setLabel('left', 'mV', **label_style)
    p.setLabel('top', f'{s}', **label_style)
    curves[s] = p.plot(pen=pen, connect='finite')  # 'finite' leaves a gap at the nan cursor instead of a line across it

sos_bp = signal.butter(4, [CARRIER_FREQ - BANDWIDTH/2, CARRIER_FREQ + BANDWIDTH/2], btype='bp', fs=Fs, output='sos')
sos_zi = {s: signal.sosfilt_zi(sos_bp) for s in SENSORS}

if LOAD_FROM_FILE:
    reader = read_from_file_realtime(LOAD_FILE_PATHS[0])

def update():
    # Sweep display: each new packet overwrites the 50 samples at the cursor and
    # nothing else moves, so previously plotted values never change.
    global write_idx
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
        for s in SENSORS:
            if FILTER:
                cleaned, sos_zi[s] = signal.sosfilt(sos_bp, data[s], zi=sos_zi[s])
                buffers[s][write_idx:write_idx + SAMPLES_PER_PACKET] = adc_to_mv(cleaned)
            else:
                buffers[s][write_idx:write_idx + SAMPLES_PER_PACKET] = adc_to_mv(data[s])
        write_idx = (write_idx + SAMPLES_PER_PACKET) % SAMPLE_MEMORY
        got_packet = True
        if LOAD_FROM_FILE:
            break  
    if got_packet:
        for s in SENSORS:
            buffers[s][write_idx:write_idx + SAMPLES_PER_PACKET] = np.nan  # blank slot marks the cursor
            curves[s].setData(buffers[s])

timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(10)
try:
    with serial_port() as ser:          # Fix 1: port always closed, even on window-close/interrupt
        if not LOAD_FROM_FILE:
            drain_to_latest(ser)
        pg.exec()
finally:
    timer.stop()  # don't leave a stale timer attached if the cell is re-run
