#!/usr/bin/env python3
"""Standalone fuzzy extractor demo for a 4-sensor photodiode array.

The script reads a CSV-like text file (such as data.txt) containing a timestamp
column plus four sensor columns. It then performs a simple enrollment/
verification flow using the same bit-extraction and Reed-Solomon logic as the
original notebook version.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

try:
    import reedsolo
except ModuleNotFoundError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "reedsolo", "--quiet"])
    import reedsolo


def bits_to_bytes(bits: Sequence[int]) -> bytes:
    """Pack a bit array into bytes, MSB-first, zero-padded to multiple of 8."""
    bits = np.asarray(bits, dtype=np.uint8)
    pad = (8 - len(bits) % 8) % 8
    bits = np.concatenate([bits, np.zeros(pad, dtype=np.uint8)])
    return bytes(
        int("".join(map(str, bits[i : i + 8])), 2)
        for i in range(0, len(bits), 8)
    )


def bytes_to_bits(b: bytes, n_bits: int | None = None) -> np.ndarray:
    """Unpack bytes back into a bit array."""
    arr = np.array([int(bit) for byte in b for bit in format(byte, "08b")], dtype=np.uint8)
    return arr[:n_bits] if n_bits is not None else arr


def extract_puf_bits(a_dict: Dict[str, float], sensors: Sequence[str], n_quant: int = 4, rank_bits: Sequence[int] | None = None) -> np.ndarray:
    """Extract rank and ratio bits from sensor amplitudes."""
    if rank_bits is not None:
        rank = np.asarray(rank_bits, dtype=np.uint8)
    else:
        rank = np.array(
            [1 if a_dict[sensors[i]] > a_dict[sensors[j]] else 0 for i in range(len(sensors)) for j in range(i + 1, len(sensors))],
            dtype=np.uint8,
        )

    geo = np.exp(np.mean([np.log(v) for v in a_dict.values()]))
    ratio: List[int] = []
    for s in sensors:
        log_r = np.log2(a_dict[s] / geo)
        clipped = np.clip((log_r + 1.0) / 2.0, 0.0, 1.0 - 1e-9)
        q = int(clipped * (2 ** n_quant))
        ratio.extend([(q >> b) & 1 for b in range(n_quant - 1, -1, -1)])

    return np.concatenate([rank, np.array(ratio, dtype=np.uint8)])


class RSFuzzyExtractor:
    """Code-offset fuzzy extractor using Reed-Solomon over GF(2^8)."""

    def __init__(self, n_puf_bytes: int, n_ecc_bytes: int = 2, key_bytes: int = 16):
        if n_puf_bytes <= n_ecc_bytes:
            raise ValueError(
                f"n_puf_bytes ({n_puf_bytes}) must exceed n_ecc_bytes ({n_ecc_bytes}). "
                f"Increase n_quant or decrease n_ecc_bytes."
            )
        self.n_puf = n_puf_bytes
        self.n_ecc = n_ecc_bytes
        self.n_data = n_puf_bytes - n_ecc_bytes
        self.n_key = key_bytes
        self.rs = reedsolo.RSCodec(n_ecc_bytes)

    @property
    def max_byte_errors(self) -> int:
        return self.n_ecc // 2

    def gen(self, puf_bits: Sequence[int]) -> Tuple[bytes, bytes]:
        """Enrollment. Returns (key, helper)."""
        puf_bytes = bits_to_bytes(puf_bits)
        if len(puf_bytes) != self.n_puf:
            raise ValueError(f"PUF gave {len(puf_bytes)} bytes, expected {self.n_puf}")
        secret = os.urandom(self.n_data)
        codeword = bytes(self.rs.encode(secret))
        helper = bytes(p ^ c for p, c in zip(puf_bytes, codeword))
        key = hashlib.sha256(secret).digest()[: self.n_key]
        return key, helper

    def rep(self, puf_bits_noisy: Sequence[int], helper: bytes, return_debug: bool = False):
        """Verification. Returns the key, and optionally the corrected PUF bytes and error positions."""
        puf_bytes = bits_to_bytes(puf_bits_noisy)
        noisy_cw = bytes(p ^ h for p, h in zip(puf_bytes, helper))
        try:
            decoded, corrected_cw, errata_pos = self.rs.decode(noisy_cw)
            key = hashlib.sha256(bytes(decoded)).digest()[: self.n_key]
            if not return_debug:
                return key
            corrected_puf_bytes = bytes(c ^ h for c, h in zip(corrected_cw, helper))
            return key, corrected_puf_bytes, errata_pos
        except reedsolo.ReedSolomonError:
            if return_debug:
                return None, None, None
            return None


def load_sensor_data(path: Path) -> Tuple[List[str], List[Dict[str, float]]]:
    """Load timestamped sensor data from a CSV-like file."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"{path} is empty or has no header row")
        all_columns = [c for c in reader.fieldnames if c is not None]
        if "timestamp" not in all_columns:
            raise ValueError(f"{path} must contain a timestamp column")
        sensor_names = [c for c in all_columns if c != "timestamp"]
        if len(sensor_names) != 4:
            raise ValueError(f"Expected 4 sensor columns, found {len(sensor_names)}")

        samples: List[Dict[str, float]] = []
        for row in reader:
            sample = {name: float(row[name]) for name in sensor_names}
            samples.append(sample)

    return sensor_names, samples


def average_window(samples: Sequence[Dict[str, float]], sensors: Sequence[str], start: int, end: int) -> Dict[str, float]:
    return {sensor: float(np.mean([samples[i][sensor] for i in range(start, end)])) for sensor in sensors}


def run_demo(data_path: Path, sensors: Sequence[str] | None = None, enroll_count: int = 20, verify_count: int = 20, window_size: int = 10, n_quant: int = 4, n_ecc_bytes: int = 2) -> None:
    sensor_names, samples = load_sensor_data(data_path)
    if sensors is None:
        sensors = tuple(sensor_names)
    else:
        sensors = tuple(sensors)
        if len(sensors) != 4:
            raise ValueError("Please provide exactly 4 sensor names")

    if len(samples) < enroll_count + verify_count:
        raise ValueError(f"Not enough samples: need at least {enroll_count + verify_count}, found {len(samples)}")

    # Enrollment uses the first window.
    enroll_window = average_window(samples, sensors, 0, enroll_count)
    rank_bits = np.array(
        [1 if enroll_window[sensors[i]] > enroll_window[sensors[j]] else 0 for i in range(len(sensors)) for j in range(i + 1, len(sensors))],
        dtype=np.uint8,
    )
    puf_bits_enroll = extract_puf_bits(enroll_window, sensors, n_quant=n_quant, rank_bits=rank_bits)
    n_bits = len(puf_bits_enroll)
    n_puf_bytes = len(bits_to_bytes(puf_bits_enroll))
    if n_ecc_bytes >= n_puf_bytes:
        n_ecc_bytes = max(1, n_puf_bytes - 1)
        print(f"Adjusted n_ecc_bytes to {n_ecc_bytes} so it stays below the {n_puf_bytes}-byte fingerprint size")

    print(f"Loaded {len(samples)} samples from {data_path}")
    print(f"Sensors: {', '.join(sensors)}")
    print(f"Enrollment window: first {enroll_count} samples")
    print(f"PUF fingerprint : {n_bits} bits → {n_puf_bytes} bytes")
    print(f"  rank bits     : {len(sensors) * (len(sensors) - 1) // 2}")
    print(f"  ratio bits    : {len(sensors) * n_quant}")
    print(f"  bit pattern   : {''.join(map(str, puf_bits_enroll))}")

    fe = RSFuzzyExtractor(n_puf_bytes=n_puf_bytes, n_ecc_bytes=n_ecc_bytes, key_bytes=16)
    print(f"RS correction capacity : {fe.max_byte_errors} bit flip(s) per read")

    key_enroll, helper = fe.gen(puf_bits_enroll)

    verify_start = enroll_count + 2
    verify_window = average_window(samples, sensors, verify_start, verify_start + verify_count)
    puf_bits_verify = extract_puf_bits(verify_window, sensors, n_quant=n_quant, rank_bits=rank_bits)
    diff = puf_bits_enroll ^ puf_bits_verify
    hd = int(diff.sum())

    print(f"Verification window: samples {verify_start}–{verify_start + verify_count - 1}")
    print(f"Enroll : {''.join(map(str, puf_bits_enroll))}")
    print(f"Verify : {''.join(map(str, puf_bits_verify))}")
    print(f"Diffs  : {''.join('X' if b else '.' for b in diff)}")
    print(f"Hamming distance : {hd}/{len(puf_bits_enroll)} bits — RS capacity {fe.max_byte_errors} bit(s)")

    key_verify, corrected_puf_bytes, errata_pos = fe.rep(puf_bits_verify, helper, return_debug=True)
    if corrected_puf_bytes is not None:
        noisy_bits = np.asarray(puf_bits_verify, dtype=np.uint8)
        corrected_bits = bytes_to_bits(corrected_puf_bytes, len(puf_bits_enroll))
        corrections = [
            (idx, int(noisy_bits[idx]), int(corrected_bits[idx]))
            for idx in range(len(noisy_bits))
            if noisy_bits[idx] != corrected_bits[idx]
        ]
        print("Corrected bit positions:")
        for idx, before, after in corrections:
            print(f"  bit {idx}: {before} -> {after}")
        print(f"RS corrected positions (codeword bytes): {list(errata_pos)}")
    else:
        print("No corrected bits were reported by the RS decoder.")

    print(f"\n{'=' * 60}")
    if key_verify is None:
        print("✗ VERIFICATION FAILED — too many bit errors for RS to correct.")
    elif key_verify != key_enroll:
        print("✗ KEY MISMATCH — RS decoded but the key does not match.")
    else:
        print("✓ VERIFICATION PASSED")
    print(f"{'=' * 60}")

    # Optional sweep over other windows to assess stability.
    max_windows = max(1, (len(samples) - verify_count - enroll_count) // window_size)
    print(f"\nSweeping {max_windows} windows of size {window_size}...")
    hds: List[int] = []
    for idx in range(max_windows):
        start = enroll_count + idx * window_size
        end = start + window_size
        if end > len(samples):
            break
        a_w = average_window(samples, sensors, start, end)
        bits_w = extract_puf_bits(a_w, sensors, n_quant=n_quant, rank_bits=rank_bits)
        diff_w = puf_bits_enroll ^ bits_w
        hd_w = int(diff_w.sum())
        hds.append(hd_w)
        key_ok, corrected_puf_bytes_w, _ = fe.rep(bits_w, helper, return_debug=True)
        if corrected_puf_bytes_w is not None:
            corrected_bits_w = bytes_to_bits(corrected_puf_bytes_w, len(puf_bits_enroll))
            corrections_w = [
                (j, int(bits_w[j]), int(corrected_bits_w[j]))
                for j in range(len(bits_w))
                if bits_w[j] != corrected_bits_w[j]
            ]
            if corrections_w:
                print(f"  window {idx + 1}: HD={hd_w}  {'PASS' if key_ok is not None else 'FAIL'}  corrected bits {corrections_w}")
                continue
        print(f"  window {idx + 1}: HD={hd_w}  {'PASS' if key_ok is not None else 'FAIL'}")

    if hds:
        print(f"HD summary: mean={np.mean(hds):.2f}, min={min(hds)}, max={max(hds)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a fuzzy-extractor test on a 4-sensor photodiode data file")
    parser.add_argument("--data", type=Path, default=Path("data.txt"), help="Path to the CSV-like input file")
    parser.add_argument("--sensors", nargs=4, default=None, help="Exactly 4 sensor column names to use")
    parser.add_argument("--enroll-count", type=int, default=20, help="Number of samples to average for enrollment")
    parser.add_argument("--verify-count", type=int, default=20, help="Number of samples to average for verification")
    parser.add_argument("--window-size", type=int, default=10, help="Window size for the stability sweep")
    parser.add_argument("--n-quant", type=int, default=4, help="Number of ratio bits per sensor")
    parser.add_argument("--n-ecc-bytes", type=int, default=2, help="Reed-Solomon error-correction bytes")
    args = parser.parse_args()

    run_demo(
        data_path=args.data,
        sensors=args.sensors,
        enroll_count=args.enroll_count,
        verify_count=args.verify_count,
        window_size=args.window_size,
        n_quant=args.n_quant,
        n_ecc_bytes=args.n_ecc_bytes,
    )


if __name__ == "__main__":
    main()
