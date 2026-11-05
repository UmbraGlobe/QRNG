r"""
## Spectrogram: Frequency Versus Time

**Purpose:** Shows how the frequency content of each channel changes over time, making intermittent tones, harmonics, and interference easier to find.

**Usage:** Set `CAPTURE_TIME` in seconds. `WINDOWSIZE` controls the time/frequency tradeoff, `OVERLAP` controls overlap between windows, and `NFFT` controls displayed frequency spacing. Leave `FILTER = False` for a full-spectrum check; enable it only when intentionally examining the carrier band.

Live mode saves `data.csv`, `summary.txt`, and `spectrograms.png` under `data/spectrograms/...`. Saved-file mode analyzes the beginning of the first file for `CAPTURE_TIME` seconds and displays the result without making a new run folder. Brighter colors mean greater power at that time and frequency.

Source: notebook cell 11 (code copied unchanged unless noted below).
"""

from qrng_common import *

### CONFIG

NFFT = 4096        
WINDOWSIZE = 512
WINDOW = 'han'
OVERLAP = 90
FILTER = False
FOLDER = "spectrograms"
CAPTURE_TIME = 30  # seconds
BANDWIDTH = 4  # Hz, for bandpass filter around carrier frequency

###

NUM_SAMPLES_TO_SAVE = int(Fs * CAPTURE_TIME)
data = read_from_file_instant(LOAD_FILE_PATHS[0]) if LOAD_FROM_FILE else read_data_realtime(NUM_SAMPLES_TO_SAVE)
data = data[:NUM_SAMPLES_TO_SAVE]  # file replay covers the same duration as live capture

win = get_window(WINDOW, WINDOWSIZE)
noverlap = int(WINDOWSIZE * OVERLAP / 100)
sos_bp = signal.butter(4, [CARRIER_FREQ - BANDWIDTH/2, CARRIER_FREQ + BANDWIDTH/2], btype='bp', fs=Fs, output='sos') 
dir = save_run_data(FOLDER, data=[data], filenames=["data.csv"], sample_count=NUM_SAMPLES_TO_SAVE, summary=f"Sampling Rate = {Fs:.2f} Hz\nCarrier Frequency = {CARRIER_FREQ:.6f} Hz\n")

with plt.style.context('dark_background'):
    fig, axes = plt.subplots(4, 1, figsize=(19, 21), sharex=False, constrained_layout=True)
    fig.suptitle(f'Channel spectrograms, Sampling Rate = {Fs:.2f} Hz', fontsize=37.5)

    for ax, s in zip(axes, SENSORS):
        if FILTER:
            spec, freqs, t_bins, im = ax.specgram(signal.sosfilt(sos_bp, adc_to_mv(data[s])), Fs=Fs, NFFT=WINDOWSIZE, window=win, noverlap=noverlap, pad_to=NFFT, cmap='magma')
        else:
            spec, freqs, t_bins, im = ax.specgram(adc_to_mv(data[s]), Fs=Fs, NFFT=WINDOWSIZE, window=win, noverlap=noverlap, pad_to=NFFT, cmap='magma')
        ax.set_ylabel(f'{s}\n(Hz)', fontsize=37.5)
        ax.tick_params(labelsize=12.5)
        ax.set_ylim(0, Fs/2)
        cbar = fig.colorbar(im, ax=ax, pad=0.01)
        cbar.set_label('Power (dB)', fontsize=37.5)
        cbar.ax.tick_params(labelsize=12.5)

    axes[-1].set_xlabel('Time (s)', fontsize=37.5)
    if dir is not None:
        plt.savefig(dir / "spectrograms.png", dpi=300)
    plt.show()

