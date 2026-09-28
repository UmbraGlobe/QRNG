
#!/usr/bin/env python3
"""
standalone_bch.py

Single-file binary BCH encoder/decoder with optional Excel/CSV batch decoding.

Change these settings for the normal workflow:

    LOAD_FILE_PATHS = [ROOT / "data" / "..." / "data_t1.csv"]
    TOTAL_BITS = 6
    BIT_CORRECTIONS = 1

TOTAL_BITS is the transmitted BCH codeword length. If it is not 2^m - 1,
the script automatically uses a shortened primitive binary BCH code.

Input data may be:
  1) raw 4-channel photodiode/ADC data with a header like:
       timestamp,A2V-16_A1,A2V-16_A2,A2V-16_A3,A2V-16_A4
     The script automatically converts each row into 6 pairwise comparison bits:
       A1>A2, A1>A3, A1>A4, A2>A3, A2>A4, A3>A4
  2) one column containing bitstrings such as 010101...
  3) TOTAL_BITS columns containing 0/1 values

For raw 4-channel data, a comparison bit is 1 when the first channel is larger
than the second channel, otherwise 0.

For repeated raw PUF fingerprints, this script performs a code-offset fuzzy
extractor analysis:

  1) Convert each raw 4-channel sample into 6 pairwise PUF bits.
  2) Enroll the modal (most common) 6-bit response.
  3) Generate BCH helper data from that enrolled response.
  4) Reconstruct every noisy response through BCH.
  5) Save a detailed CSV and a PNG showing:
       - per-bit stability
       - intra-Hamming-distance distribution
       - raw exact-match rate
       - BCH recovery rate
       - corrected/failing reads
"""

from __future__ import annotations

import csv
import os
import random
import textwrap
from collections import Counter
from array import array
from pathlib import Path
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple


# ============================================================
# USER SETTINGS
# ============================================================

ROOT = Path(__file__).resolve().parent

LOAD_FILE_PATHS = [
    ROOT / "data/wavelength_stream/08_04_2026/green/data.csv",
]

#  d_min >= 2t + 1 ( t is the bit corrections)
# Raw 4-channel input settings.
TOTAL_BITS = 6              # Six pairwise comparison bits from four sensors
BIT_CORRECTIONS = 1         # Maximum guaranteed bit errors corrected, t

# Raw 4-channel input settings.
# These match the header in your data.csv file.
SENSOR_COLUMNS = [
    "A2V-16_A1",
    "A2V-16_A2",
    "A2V-16_A3",
    "A2V-16_A4",
]

# Pair order used to make the 6-bit PUF response:
# bit0=A1>A2, bit1=A1>A3, bit2=A1>A4,
# bit3=A2>A3, bit4=A2>A4, bit5=A3>A4
PAIR_INDICES = [
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 2),
    (1, 3),
    (2, 3),
]

# Optional settings
SHEET_NAME = None           # None = active Excel sheet
OUTPUT_DIR = None           # None = save results beside each input file
RUN_SELF_TEST = True
SELF_TEST_TRIALS = 200
RANDOM_SEED = 0

# Visualization / analysis settings
SAVE_SUMMARY_PNG = True

# Match the publication-style visualization theme used by the analysis figures.
SUMMARY_WIDTH_IN = 7.16
SUMMARY_DPI = 150
SUMMARY_SAVE_DPI = 600
SUMMARY_SERIF = ["Times New Roman", "Nimbus Roman", "DejaVu Serif"]

# Okabe-Ito qualitative palette -- same fixed order as the analysis visualization code.
SUMMARY_SERIES = [
    "#0072B2", "#D55E00", "#009E73", "#E69F00",
    "#CC79A7", "#56B4E9", "#F0E442", "#000000",
]
SUMMARY_STATUS = {
    "good": "#000000",
    "warning": "#D55E00",
    "serious": "#D55E00",
    "critical": "#B2182B",
}
SUMMARY_INK = {
    "ink": "#000000",
    "ink2": "#333333",
    "muted": "#606060",
    "rule": "#000000",
    "grid": "#D9D9D9",
    "surface": "#FFFFFF",
}

SUMMARY_RC = {
    "font.family": "serif",
    "font.serif": SUMMARY_SERIF,
    "mathtext.fontset": "stix",
    "font.size": 14,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "axes.linewidth": 0.6,
    "axes.edgecolor": "#000000",
    "axes.labelcolor": "#000000",
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.color": "#000000",
    "ytick.color": "#000000",
    "lines.linewidth": 1.0,
    "lines.markersize": 3.5,
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "legend.columnspacing": 1.2,
    "legend.labelspacing": 0.35,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}

# Exact rolling-average window retained from the BCH analysis.
ROLLING_MEAN_READS = 601

# Display downsampling only; statistics and rolling mean still use every read.
MAX_SCATTER_POINTS = 150_000

PAIR_LABELS_SLASH = [
    "A1/A2",
    "A1/A3",
    "A1/A4",
    "A2/A3",
    "A2/A4",
    "A3/A4",
]


# ============================================================
# GF(2^m) / BCH IMPLEMENTATION
# ============================================================

# Primitive polynomials represented as integers, including x^m and constant 1.
# These are standard primitive choices for GF(2^m).
PRIMITIVE_POLYNOMIALS = {
    2:  0b111,                 # x^2 + x + 1
    3:  0b1011,                # x^3 + x + 1
    4:  0b10011,               # x^4 + x + 1
    5:  0b100101,              # x^5 + x^2 + 1
    6:  0b1000011,             # x^6 + x + 1
    7:  0b10000011,            # x^7 + x + 1
    8:  0b100011101,           # x^8 + x^4 + x^3 + x^2 + 1
    9:  0b1000010001,          # x^9 + x^4 + 1
    10: 0b10000001001,         # x^10 + x^3 + 1
    11: 0b100000000101,        # x^11 + x^2 + 1
    12: 0b1000010100111,       # x^12 + x^6 + x^4 + x + 1
    13: 0b10000000011011,      # x^13 + x^4 + x^3 + x + 1
    14: 0b100010001000011,     # x^14 + x^10 + x^6 + x + 1
    15: 0b1000000000000011,    # x^15 + x + 1
}


def _poly_degree(poly: int) -> int:
    return poly.bit_length() - 1


def _gf2_mod(dividend: int, divisor: int) -> int:
    """Polynomial remainder over GF(2), both polynomials stored as bitmasks."""
    ddeg = _poly_degree(divisor)
    while dividend and _poly_degree(dividend) >= ddeg:
        dividend ^= divisor << (_poly_degree(dividend) - ddeg)
    return dividend


class GF2m:
    def __init__(self, m: int, primitive_polynomial: int):
        self.m = m
        self.n = (1 << m) - 1
        self.primitive_polynomial = primitive_polynomial

        if _poly_degree(primitive_polynomial) != m:
            raise ValueError("Primitive polynomial degree does not match m.")

        self.exp = [0] * (2 * self.n)
        self.log = [-1] * (1 << m)

        x = 1
        for i in range(self.n):
            if self.log[x] != -1:
                raise ValueError(
                    f"Polynomial 0x{primitive_polynomial:x} is not primitive for m={m}."
                )
            self.exp[i] = x
            self.log[x] = i

            x <<= 1
            if x & (1 << m):
                x ^= primitive_polynomial
            x &= self.n

        if x != 1:
            raise ValueError(
                f"Polynomial 0x{primitive_polynomial:x} is not primitive for m={m}."
            )

        for i in range(self.n, 2 * self.n):
            self.exp[i] = self.exp[i - self.n]

    def add(self, a: int, b: int) -> int:
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def div(self, a: int, b: int) -> int:
        if b == 0:
            raise ZeroDivisionError("GF division by zero")
        if a == 0:
            return 0
        return self.exp[(self.log[a] - self.log[b]) % self.n]

    def pow_alpha(self, exponent: int) -> int:
        return self.exp[exponent % self.n]


def _poly_mul_gf(a: Sequence[int], b: Sequence[int], gf: GF2m) -> List[int]:
    """Multiply low-order-first polynomials whose coefficients are GF elements."""
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj:
                out[i + j] ^= gf.mul(ai, bj)
    return out


def _cyclotomic_coset(start: int, n: int) -> Tuple[int, ...]:
    seen = []
    x = start % n
    while x not in seen:
        seen.append(x)
        x = (x * 2) % n
    return tuple(seen)


def _generator_polynomial(gf: GF2m, t: int) -> Tuple[int, List[Tuple[int, ...]]]:
    """
    Build the narrow-sense primitive binary BCH generator polynomial having
    alpha^1 ... alpha^(2t) as roots.
    """
    n = gf.n
    used = set()
    cosets = []

    for exponent in range(1, 2 * t + 1):
        coset = _cyclotomic_coset(exponent, n)
        key = min(coset)
        if key not in used:
            used.add(key)
            cosets.append(coset)

    g = [1]  # coefficients low-order first, initially polynomial 1

    for coset in cosets:
        minimal = [1]
        for exponent in coset:
            # x + alpha^exponent  (same as x - alpha^exponent in characteristic 2)
            minimal = _poly_mul_gf(minimal, [gf.pow_alpha(exponent), 1], gf)

        # A binary minimal polynomial must have coefficients only 0 or 1.
        if any(c not in (0, 1) for c in minimal):
            raise RuntimeError(
                f"Internal error: minimal polynomial for coset {coset} "
                "did not reduce to GF(2)."
            )

        # g is binary so multiplying in the extension field is fine.
        g = _poly_mul_gf(g, minimal, gf)

    if any(c not in (0, 1) for c in g):
        raise RuntimeError("Internal error: BCH generator is not binary.")

    generator = 0
    for degree, coeff in enumerate(g):
        if coeff:
            generator |= 1 << degree

    return generator, cosets


@dataclass
class DecodeResult:
    success: bool
    corrected: int
    message: int
    error_positions: List[int]
    reason: str = ""


class BCH:
    """
    Primitive narrow-sense binary BCH code, shortened from n_full = 2^m - 1
    when requested n is smaller.
    """

    def __init__(self, n: int, t: int):
        if n < 3:
            raise ValueError("TOTAL_BITS must be at least 3.")
        if t < 1:
            raise ValueError("BIT_CORRECTIONS must be at least 1.")

        self.n = int(n)
        self.t = int(t)

        m = 2
        while (1 << m) - 1 < self.n:
            m += 1

        if m not in PRIMITIVE_POLYNOMIALS:
            raise ValueError(
                f"TOTAL_BITS={self.n} needs m={m}, but this script supports "
                f"m={min(PRIMITIVE_POLYNOMIALS)}..{max(PRIMITIVE_POLYNOMIALS)}."
            )

        self.m = m
        self.n_full = (1 << m) - 1
        self.shorten = self.n_full - self.n
        self.gf = GF2m(m, PRIMITIVE_POLYNOMIALS[m])

        self.generator, self.cosets = _generator_polynomial(self.gf, self.t)
        self.parity_bits = _poly_degree(self.generator)
        self.k_full = self.n_full - self.parity_bits
        self.k = self.n - self.parity_bits

        if self.k <= 0:
            raise ValueError(
                f"Requested n={self.n}, t={self.t} leaves no message bits. "
                f"The BCH generator needs {self.parity_bits} parity bits. "
                "Increase TOTAL_BITS or reduce BIT_CORRECTIONS."
            )

    def encode(self, message: int) -> int:
        if message < 0 or message >= (1 << self.k):
            raise ValueError(f"Message must fit in k={self.k} bits.")

        shifted = message << self.parity_bits
        remainder = _gf2_mod(shifted, self.generator)
        codeword = shifted ^ remainder

        if codeword >= (1 << self.n):
            raise RuntimeError(
                "Shortening error: generated codeword did not fit requested length."
            )
        return codeword

    def syndromes(self, received: int) -> List[int]:
        if received < 0 or received >= (1 << self.n_full):
            raise ValueError("Received word is too large.")

        set_positions = [p for p in range(self.n_full) if (received >> p) & 1]
        synd = []
        for j in range(1, 2 * self.t + 1):
            s = 0
            for p in set_positions:
                s ^= self.gf.pow_alpha(j * p)
            synd.append(s)
        return synd

    def _berlekamp_massey(self, synd: Sequence[int]) -> Tuple[List[int], int]:
        N = len(synd)
        C = [0] * (N + 1)
        B = [0] * (N + 1)
        C[0] = 1
        B[0] = 1

        L = 0
        shift = 1
        b = 1

        for n_idx in range(N):
            discrepancy = synd[n_idx]
            for i in range(1, L + 1):
                if C[i] and synd[n_idx - i]:
                    discrepancy ^= self.gf.mul(C[i], synd[n_idx - i])

            if discrepancy == 0:
                shift += 1
                continue

            old_C = C[:]
            scale = self.gf.div(discrepancy, b)

            for j in range(0, N + 1 - shift):
                if B[j]:
                    C[j + shift] ^= self.gf.mul(scale, B[j])

            if 2 * L <= n_idx:
                L = n_idx + 1 - L
                B = old_C
                b = discrepancy
                shift = 1
            else:
                shift += 1

        return C[: L + 1], L

    def _chien_search(self, locator: Sequence[int], degree: int) -> List[int]:
        positions = []

        for p in range(self.n_full):
            # Error at coefficient x^p => locator root at alpha^(-p)
            x = self.gf.pow_alpha(-p)
            value = 0
            x_power = 1
            for coeff in locator:
                if coeff:
                    value ^= self.gf.mul(coeff, x_power)
                x_power = self.gf.mul(x_power, x)

            if value == 0:
                positions.append(p)

        if len(positions) != degree:
            return []

        return positions

    def decode(self, received: int) -> DecodeResult:
        if received < 0 or received >= (1 << self.n):
            return DecodeResult(
                False, received, 0, [], f"Received word must be exactly {self.n} bits or less."
            )

        # Shortened codewords are restored to the parent code by prepending
        # self.shorten zeros. Integer representation already does this.
        full_received = received

        synd = self.syndromes(full_received)
        if not any(synd):
            return DecodeResult(
                True,
                received,
                received >> self.parity_bits,
                [],
                "No errors detected.",
            )

        locator, degree = self._berlekamp_massey(synd)

        if degree == 0 or degree > self.t:
            return DecodeResult(
                False,
                received,
                received >> self.parity_bits,
                [],
                f"Decoder estimated {degree} errors; correction limit is t={self.t}.",
            )

        positions = self._chien_search(locator, degree)
        if len(positions) != degree:
            return DecodeResult(
                False,
                received,
                received >> self.parity_bits,
                [],
                "Could not locate a consistent set of error positions.",
            )

        # Positions removed during shortening were never transmitted.
        if any(p >= self.n for p in positions):
            return DecodeResult(
                False,
                received,
                received >> self.parity_bits,
                positions,
                "Decoder located an error in a shortened-away bit position.",
            )

        corrected = full_received
        for p in positions:
            corrected ^= 1 << p

        if any(self.syndromes(corrected)):
            return DecodeResult(
                False,
                received,
                received >> self.parity_bits,
                positions,
                "Syndromes were nonzero after attempted correction.",
            )

        corrected_short = corrected & ((1 << self.n) - 1)
        message = corrected_short >> self.parity_bits

        return DecodeResult(
            True,
            corrected_short,
            message,
            sorted(positions),
            f"Corrected {len(positions)} bit error(s).",
        )

    def summary(self) -> str:
        return (
            f"BCH parameters\n"
            f"  requested codeword bits n = {self.n}\n"
            f"  guaranteed corrections t  = {self.t}\n"
            f"  field GF(2^{self.m})\n"
            f"  parent BCH length          = {self.n_full}\n"
            f"  shortened bits removed     = {self.shorten}\n"
            f"  message bits k             = {self.k}\n"
            f"  parity bits n-k            = {self.parity_bits}\n"
            f"  generator polynomial       = 0x{self.generator:x}"
        )


# ============================================================
# FILE LOADING
# ============================================================

def _normalize_bitstring(value) -> Optional[str]:
    if value is None:
        return None

    if isinstance(value, bool):
        return "1" if value else "0"

    if isinstance(value, int) and value in (0, 1):
        return str(value)

    if isinstance(value, float) and value in (0.0, 1.0):
        return str(int(value))

    text = str(value).strip().replace(" ", "").replace("_", "")
    if text.startswith("'"):
        text = text[1:]

    if text and all(ch in "01" for ch in text):
        return text
    return None


def _extract_word_from_row(row: Sequence, n: int) -> Optional[str]:
    """
    Fallback parser for already-generated bit data:
      - first n individual 0/1 cells found in the row, OR
      - one cell containing an n-character bitstring.
    """
    individual = []
    bitstrings = []

    for value in row:
        bits = _normalize_bitstring(value)
        if bits is None:
            continue

        if len(bits) == 1:
            individual.append(bits)
        elif len(bits) == n:
            bitstrings.append(bits)

    if len(individual) >= n:
        return "".join(individual[:n])

    if bitstrings:
        return bitstrings[0]

    return None


def _pairwise_bits(values: Sequence[float]) -> str:
    """
    Convert four sensor values into six pairwise comparison bits.

    Order:
      0: A1 > A2
      1: A1 > A3
      2: A1 > A4
      3: A2 > A3
      4: A2 > A4
      5: A3 > A4

    Ties produce 0.
    """
    if len(values) != 4:
        raise ValueError("Exactly four sensor values are required.")

    return "".join(
        "1" if values[i] > values[j] else "0"
        for i, j in PAIR_INDICES
    )


def _find_sensor_indices(header: Sequence) -> Optional[List[int]]:
    """
    Return the four sensor-column indices when the expected raw-data header
    is present. Matching is case-insensitive and ignores surrounding spaces.
    """
    normalized = {
        str(value).strip().lower(): idx
        for idx, value in enumerate(header)
        if value is not None
    }

    wanted = [name.strip().lower() for name in SENSOR_COLUMNS]
    if all(name in normalized for name in wanted):
        return [normalized[name] for name in wanted]

    return None


def _word_from_analog_row(row: Sequence, sensor_indices: Sequence[int]) -> Optional[str]:
    """Read four raw channel values and create the six comparison bits."""
    try:
        values = [float(row[idx]) for idx in sensor_indices]
    except (TypeError, ValueError, IndexError):
        return None

    return _pairwise_bits(values)


def iter_words(path, n: int, sheet_name=None):
    """
    Yield n-bit words one row at a time.

    Auto-detection:
      1) If the file contains SENSOR_COLUMNS, interpret it as raw 4-channel
         analog/ADC data and generate six pairwise comparison bits per row.
      2) Otherwise fall back to parsing existing bitstrings / 0-1 columns.

    Streaming keeps memory usage low for million-row CSV files.
    """
    path = Path(path)
    ext = path.suffix.lower()

    if ext == ".csv":
        with path.open("r", newline="", encoding="utf-8-sig") as f:
            reader = csv.reader(f)

            try:
                first_row = next(reader)
            except StopIteration:
                return

            sensor_indices = _find_sensor_indices(first_row)

            if sensor_indices is not None:
                if n != 6:
                    raise ValueError(
                        "Raw 4-channel mode generates exactly 6 pairwise bits, "
                        f"but TOTAL_BITS={n}. Set TOTAL_BITS = 6."
                    )

                for row in reader:
                    bits = _word_from_analog_row(row, sensor_indices)
                    if bits is not None:
                        yield bits
                return

            # No raw-data header detected: treat first row and the rest as
            # already-generated bit data.
            bits = _extract_word_from_row(first_row, n)
            if bits is not None:
                yield bits

            for row in reader:
                bits = _extract_word_from_row(row, n)
                if bits is not None:
                    yield bits
        return

    if ext in (".xlsx", ".xlsm"):
        try:
            from openpyxl import load_workbook
        except ImportError as exc:
            raise SystemExit(
                "Reading .xlsx files requires openpyxl. Install it with:\n"
                "    pip install openpyxl"
            ) from exc

        wb = load_workbook(path, read_only=True, data_only=True)
        try:
            ws = wb[sheet_name] if sheet_name else wb.active
            rows = ws.iter_rows(values_only=True)

            try:
                first_row = next(rows)
            except StopIteration:
                return

            sensor_indices = _find_sensor_indices(first_row)

            if sensor_indices is not None:
                if n != 6:
                    raise ValueError(
                        "Raw 4-channel mode generates exactly 6 pairwise bits, "
                        f"but TOTAL_BITS={n}. Set TOTAL_BITS = 6."
                    )

                for row in rows:
                    bits = _word_from_analog_row(row, sensor_indices)
                    if bits is not None:
                        yield bits
                return

            bits = _extract_word_from_row(first_row, n)
            if bits is not None:
                yield bits

            for row in rows:
                bits = _extract_word_from_row(row, n)
                if bits is not None:
                    yield bits
        finally:
            wb.close()
        return

    raise ValueError("Input file must end in .csv, .xlsx, or .xlsm")


# ============================================================
# TEST / BATCH RUNNER
# ============================================================

def self_test(code: BCH, trials: int = 200, seed: int = 0) -> Tuple[int, int]:
    rng = random.Random(seed)
    passed = 0

    for _ in range(trials):
        msg = rng.getrandbits(code.k)
        original = code.encode(msg)

        error_count = rng.randint(0, code.t)
        error_positions = rng.sample(range(code.n), error_count)

        received = original
        for pos in error_positions:
            received ^= 1 << pos

        result = code.decode(received)

        if (
            result.success
            and result.corrected == original
            and result.message == msg
        ):
            passed += 1

    return passed, trials


def output_path_for(input_path, suffix: str) -> Path:
    input_path = Path(input_path)

    if OUTPUT_DIR is None:
        return input_path.with_name(input_path.stem + suffix)

    output_dir = Path(OUTPUT_DIR)
    if not output_dir.is_absolute():
        output_dir = ROOT / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / (input_path.stem + suffix)


def default_output_path(input_path) -> Path:
    return output_path_for(input_path, "_bch_results.csv")


def _hamming_distance(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def _int_to_bits(value: int, n: int) -> str:
    return format(value, f"0{n}b")


def _number_word(value: int, title=False) -> str:
    words = {
        0: "zero",
        1: "one",
        2: "two",
        3: "three",
        4: "four",
        5: "five",
        6: "six",
        7: "seven",
        8: "eight",
        9: "nine",
        10: "ten",
    }
    out = words.get(value, str(value))
    return out.title() if title else out


def _rolling_mean_centered(values, window: int):
    """
    O(N) centered rolling mean with NaNs at the two incomplete edges.
    This is fast enough for ~1,000,000 fingerprint reads.
    """
    import numpy as np

    x = np.asarray(values, dtype=float)
    n = len(x)

    if n == 0:
        return x

    window = max(1, int(window))

    if window > n:
        window = n

    # Keep the requested window exactly as specified.
    # For an even centered window (e.g. 6000), the alignment differs by one
    # sample between the left and right edges, but every valid mean still
    # contains exactly `window` reads.
    if window <= 1:
        return x.copy()

    csum = np.empty(n + 1, dtype=float)
    csum[0] = 0.0
    np.cumsum(x, out=csum[1:])

    valid = (csum[window:] - csum[:-window]) / window

    left = window // 2
    right = n - len(valid) - left

    return np.concatenate(
        (
            np.full(left, np.nan),
            valid,
            np.full(right, np.nan),
        )
    )


def _summary_axes_style(ax, T, grid=False, gridaxis="y"):
    """Apply the same boxed-axis styling as the publication visualization code."""
    ylabel = ax.get_ylabel()
    if len(ylabel) > 24 and "\n" not in ylabel:
        ax.set_ylabel(textwrap.fill(ylabel, width=22, break_long_words=False,
                                     break_on_hyphens=False))
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)
        spine.set_color(T["rule"])
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    if grid:
        ax.grid(axis=gridaxis, color=T["grid"], lw=0.4, ls=":")
        ax.set_axisbelow(True)


def _summary_legend(ax, T, n=2, frame=True, **kwargs):
    """Use the same compact journal-style legend as the other visualization code."""
    kwargs.setdefault("loc", "best")
    legend = ax.legend(
        frameon=frame,
        facecolor="white",
        edgecolor=T["rule"] if frame else "none",
        framealpha=1.0,
        fancybox=False,
        borderpad=0.35,
        labelspacing=0.3,
        labelcolor=T["ink"],
        fontsize=11,
        ncols=min(n, 2),
        handlelength=1.4,
        columnspacing=1.0,
        **kwargs,
    )
    if frame:
        legend.get_frame().set_linewidth(0.4)
    return legend


def _summary_panel_title(ax, letter, title, T):
    """Give each panel letter a little more emphasis than its title."""
    for label, x_offset, font_size in ((f"({letter})", 0, 16), (title, 29, 14)):
        ax.annotate(
            label,
            xy=(0, 1),
            xycoords="axes fraction",
            xytext=(x_offset, 3.5),
            textcoords="offset points",
            color=T["ink"],
            fontsize=font_size,
            ha="left",
            va="bottom",
        )


def _make_summary_png(
    file_path: Path,
    code: BCH,
    enrolled: str,
    helper_data: str,
    total: int,
    raw_exact: int,
    recovered: int,
    mismatch_counts: List[int],
    corrected_bit_counts: List[int],
    hd_counts: List[int],
    mean_hd_bits: float,
    hd_values,
) -> Path:
    """
    Graph-only BCH recovery visualization.

    The summary/statistics strip and caption are intentionally omitted so the
    exported figure matches the clean publication-style layout used by the
    reference figure: only plotted data panels are shown.
    """
    try:
        import numpy as np
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "Creating the PNG requires numpy and matplotlib.\n"
            "Install them with:\n"
            "    pip install numpy matplotlib"
        ) from exc

    png_path = output_path_for(file_path, "_bch_recovery_plots.png")

    # -----------------------------------------------------------------
    # Plot data
    # -----------------------------------------------------------------
    raw_flip_pct = np.asarray(mismatch_counts, dtype=float) * 100.0 / total
    fixed_flip_pct = np.asarray(corrected_bit_counts, dtype=float) * 100.0 / total
    hd_dist_pct = np.asarray(hd_counts, dtype=float) * 100.0 / total

    hd_values_np = np.asarray(hd_values, dtype=np.uint8)
    per_read_hd_pct = hd_values_np.astype(float) * (100.0 / code.n)
    rolling = _rolling_mean_centered(per_read_hd_pct, ROLLING_MEAN_READS)
    reads = np.arange(1, total + 1)

    T = dict(SUMMARY_INK)
    T["series"] = SUMMARY_SERIES
    BLUE = T["series"][0]
    ORANGE = T["series"][1]
    GREEN = T["series"][2]
    RED = SUMMARY_STATUS["critical"]
    GRAY = "#BDBDBD"

    # Full-width vertically stacked panels, matching the reference figure style.
    fig_height = 8.5

    with plt.rc_context(SUMMARY_RC):
        fig = plt.figure(
            figsize=(SUMMARY_WIDTH_IN, fig_height),
            dpi=SUMMARY_DPI,
        )
        gs = fig.add_gridspec(
            3,
            1,
            height_ratios=[1.55, 1.20, 1.20],
            left=0.095,
            right=0.985,
            top=0.975,
            bottom=0.065,
            hspace=0.52,
        )

        # =============================================================
        # (a) Intra-Hamming Distance Across Reads
        # =============================================================
        ax_a = fig.add_subplot(gs[0, 0])
        _summary_panel_title(ax_a, "a", "Intra-Hamming Distance Across Reads", T)

        stride = (
            int(np.ceil(total / MAX_SCATTER_POINTS))
            if total > MAX_SCATTER_POINTS
            else 1
        )
        ax_a.scatter(
            reads[::stride],
            per_read_hd_pct[::stride],
            s=1.0,
            c=GRAY,
            alpha=0.28,
            linewidths=0,
            rasterized=True,
            label="Per-read intra-HD",
            zorder=1,
        )
        ax_a.plot(
            reads,
            rolling,
            color=BLUE,
            lw=1.0,
            label=f"Rolling mean ({ROLLING_MEAN_READS} reads)",
            zorder=3,
        )

        correction_radius_pct = 100.0 * code.t / code.n
        radius_label = (
            "One-bit BCH correction radius"
            if code.t == 1
            else f"{code.t}-bit BCH correction radius"
        )
        ax_a.axhline(
            correction_radius_pct,
            color=ORANGE,
            lw=0.7,
            ls=(0, (4, 2.5)),
            label=radius_label,
            zorder=2,
        )

        max_seen_pct = float(np.max(per_read_hd_pct)) if total else 0.0
        y_top = max(50.0, np.ceil(max_seen_pct / 10.0) * 10.0)
        y_top = min(100.0, y_top)
        if y_top <= 0:
            y_top = 50.0

        ax_a.set_ylim(0, y_top + 2.5)
        ax_a.set_xlim(1, min(total, 100_000))
        ax_a.set_ylabel("Distance from Enrolled Response (%)")
        ax_a.set_xlabel("Fingerprint Read")
        _summary_axes_style(ax_a, T, grid=True, gridaxis="y")
        _summary_legend(ax_a, T, n=2, frame=True, loc="upper right")

        # =============================================================
        # (b) Bit Positions Fixed by BCH
        # =============================================================
        ax_b = fig.add_subplot(gs[1, 0])
        _summary_panel_title(ax_b, "b", "Bit Positions Fixed by BCH", T)

        x = np.arange(code.n)
        width = 0.37

        ax_b.bar(
            x - width / 2,
            raw_flip_pct,
            width=width,
            facecolor=BLUE,
            edgecolor=T["rule"],
            lw=0.5,
            label="Raw flips",
            zorder=3,
        )
        fixed_bars = ax_b.bar(
            x + width / 2,
            fixed_flip_pct,
            width=width,
            facecolor=GREEN,
            edgecolor=T["rule"],
            lw=0.5,
            label="Fixed by BCH",
            zorder=3,
        )

        pair_labels = []
        for i in range(code.n):
            pair = PAIR_LABELS_SLASH[i] if i < len(PAIR_LABELS_SLASH) else ""
            pair_labels.append(f"B{i + 1}\n{pair}")

        ax_b.set_xticks(x)
        ax_b.set_xticklabels(pair_labels, linespacing=1.1)
        ax_b.tick_params(axis="x", which="minor", bottom=False, top=False)
        ax_b.set_ylabel("Reads (%)")

        max_bit_pct = max(
            float(np.max(raw_flip_pct)) if code.n else 0.0,
            float(np.max(fixed_flip_pct)) if code.n else 0.0,
        )
        ax_b.set_ylim(0, max(5.0, max_bit_pct * 1.20))
        _summary_axes_style(ax_b, T, grid=True, gridaxis="y")
        _summary_legend(ax_b, T, n=2, frame=True, loc="upper right")

        for bar, count in zip(fixed_bars, corrected_bit_counts):
            if count <= 0:
                continue
            ax_b.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(0.12, max_bit_pct * 0.012),
                f"{count:,}",
                ha="center",
                va="bottom",
                fontsize=10,
                color=T["ink2"],
            )

        # =============================================================
        # (c) Intra-Hamming Distance Distribution
        # =============================================================
        ax_c = fig.add_subplot(gs[2, 0])
        _summary_panel_title(ax_c, "c", "Intra-Hamming Distance Distribution", T)

        hd_x = np.arange(code.n + 1)
        hd_colors = [GREEN if distance <= code.t else RED for distance in hd_x]
        hd_bars = ax_c.bar(
            hd_x,
            hd_dist_pct,
            width=0.62,
            facecolor=hd_colors,
            edgecolor=T["rule"],
            lw=0.5,
            zorder=3,
        )

        ax_c.set_xticks(hd_x)
        ax_c.set_xlabel(f"Differing Bits out of {code.n}")
        ax_c.set_ylabel("Reads (%)")

        max_hd_pct = float(np.max(hd_dist_pct)) if len(hd_dist_pct) else 0.0
        hd_y_top = max(10.0, max_hd_pct * 1.14)
        ax_c.set_ylim(0, hd_y_top)
        _summary_axes_style(ax_c, T, grid=True, gridaxis="y")

        limit_x = code.t + 0.5
        ax_c.axvline(
            limit_x,
            color=ORANGE,
            lw=0.7,
            ls=(0, (4, 2.5)),
            zorder=4,
        )
        ax_c.annotate(
            f"BCH $t$ = {code.t} limit",
            xy=(limit_x, hd_y_top * 0.82),
            xytext=(3, 0),
            textcoords="offset points",
            color=T["ink2"],
            fontsize=10,
            ha="left",
            va="center",
        )

        for bar, pct, count in zip(hd_bars, hd_dist_pct, hd_counts):
            if count <= 0:
                continue
            ax_c.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + hd_y_top * 0.012,
                f"{pct:.2f}%\n({count:,})",
                ha="center",
                va="bottom",
                fontsize=10,
                color=T["ink2"],
                linespacing=1.05,
            )

        fig.savefig(
            png_path,
            dpi=SUMMARY_SAVE_DPI,
            bbox_inches="tight",
            pad_inches=0.02,
        )
        plt.close(fig)

    return png_path


def process_file(file_path, code: BCH) -> None:
    file_path = Path(file_path)

    if not file_path.is_absolute():
        file_path = ROOT / file_path

    if not file_path.is_file():
        print(f"\n[SKIP] File not found: {file_path}")
        return

    output_path = default_output_path(file_path)

    print(f"\nReading: {file_path}")
    print("  Raw 4-channel files are converted to 6 bits using all unique sensor pairs.")

    # ========================================================
    # PASS 1: determine the enrolled/modal PUF response
    # ========================================================
    distribution = Counter()
    total = 0

    for bitstring in iter_words(file_path, code.n, SHEET_NAME):
        distribution[bitstring] += 1
        total += 1

    if total == 0:
        print(
            f"\n[SKIP] No {code.n}-bit rows were found in:\n{file_path}\n"
            "Expected either:\n"
            "  - raw columns: timestamp + A2V-16_A1 ... A2V-16_A4, or\n"
            "  - one bitstring column, or\n"
            "  - one 0/1 column per bit."
        )
        return

    # Highest-frequency response is the enrollment fingerprint.
    # Lexicographical ordering gives deterministic behavior if two responses tie.
    enrolled, enrolled_count = sorted(
        distribution.items(),
        key=lambda item: (-item[1], item[0]),
    )[0]

    # ========================================================
    # CODE-OFFSET FUZZY EXTRACTOR ENROLLMENT
    # ========================================================
    #
    # Pick a deterministic valid BCH codeword and store:
    #
    #       helper = enrolled_PUF XOR BCH_codeword
    #
    # For a noisy read:
    #
    #       candidate = noisy_PUF XOR helper
    #                 = BCH_codeword XOR PUF_error
    #
    # If PUF_error contains <= t flipped bits, BCH can recover the
    # original codeword and therefore reconstruct the enrolled PUF.
    #
    rng = random.Random(RANDOM_SEED)
    enrollment_message = rng.getrandbits(code.k)
    enrollment_codeword_int = code.encode(enrollment_message)
    enrollment_codeword = _int_to_bits(enrollment_codeword_int, code.n)

    enrolled_int = int(enrolled, 2)
    helper_int = enrolled_int ^ enrollment_codeword_int
    helper_data = _int_to_bits(helper_int, code.n)

    print(f"  enrolled/modal response: {enrolled}")
    print(f"  modal count:              {enrolled_count:,}/{total:,}")
    print(f"  enrollment BCH codeword: {enrollment_codeword}")
    print(f"  helper data:             {helper_data}")

    # ========================================================
    # PASS 2: reconstruct every read and collect statistics
    # ========================================================
    header = [
        "row",
        "received_puf",
        "hamming_distance_to_enrolled",
        "raw_exact_match",
        "helper_data",
        "bch_candidate",
        "decoder_success",
        "recovered_enrolled_response",
        "num_corrected",
        "corrected_puf_bit_positions",
        "decoded_codeword",
        "reason",
    ]

    raw_exact = 0
    recovered = 0
    mismatch_counts = [0] * code.n
    corrected_bit_counts = [0] * code.n
    hd_counts = [0] * (code.n + 1)
    hd_sum = 0

    # One unsigned byte per fingerprint read: ~1 MB for one million rows.
    # Needed only for panel (b), the per-read / rolling intra-HD visualization.
    hd_values = array("B")

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for row_index, bitstring in enumerate(
            iter_words(file_path, code.n, SHEET_NAME),
            start=1,
        ):
            received_int = int(bitstring, 2)

            hd = _hamming_distance(bitstring, enrolled)
            hd_counts[hd] += 1
            hd_sum += hd
            hd_values.append(hd)

            raw_match = bitstring == enrolled
            if raw_match:
                raw_exact += 1

            mismatch_positions = [
                i
                for i, (received_bit, enrolled_bit)
                in enumerate(zip(bitstring, enrolled))
                if received_bit != enrolled_bit
            ]

            for i in mismatch_positions:
                mismatch_counts[i] += 1

            candidate_int = received_int ^ helper_int
            result = code.decode(candidate_int)

            # Do not count a miscorrection to another valid BCH codeword.
            # Recovery is successful only if the decoder returns the exact
            # codeword used during enrollment.
            recovery_success = (
                result.success
                and result.corrected == enrollment_codeword_int
            )

            if recovery_success:
                recovered += 1

                for i in mismatch_positions:
                    corrected_bit_counts[i] += 1

            writer.writerow([
                row_index,
                bitstring,
                hd,
                raw_match,
                helper_data,
                _int_to_bits(candidate_int, code.n),
                result.success,
                recovery_success,
                len(result.error_positions) if recovery_success else "",
                " ".join(f"B{i + 1}" for i in mismatch_positions)
                if recovery_success
                else "",
                _int_to_bits(result.corrected, code.n),
                result.reason,
            ])

    mean_hd_bits = hd_sum / total
    raw_pct = 100.0 * raw_exact / total
    recovery_pct = 100.0 * recovered / total

    print(f"\nProcessed: {file_path}")
    print(f"  rows:                         {total:,}")
    print(f"  enrolled/modal response:      {enrolled}")
    print(f"  raw exact matches:            {raw_exact:,}/{total:,} ({raw_pct:.4f}%)")
    print(
        f"  BCH recovered (t={code.t}):          "
        f"{recovered:,}/{total:,} ({recovery_pct:.4f}%)"
    )
    print(
        f"  mean intra-HD:                "
        f"{mean_hd_bits:.4f}/{code.n} "
        f"({100.0 * mean_hd_bits / code.n:.4f}%)"
    )
    print(f"  results CSV:                  {output_path}")

    if SAVE_SUMMARY_PNG:
        png_path = _make_summary_png(
            file_path=file_path,
            code=code,
            enrolled=enrolled,
            helper_data=helper_data,
            total=total,
            raw_exact=raw_exact,
            recovered=recovered,
            mismatch_counts=mismatch_counts,
            corrected_bit_counts=corrected_bit_counts,
            hd_counts=hd_counts,
            mean_hd_bits=mean_hd_bits,
            hd_values=hd_values,
        )
        print(f"  plots PNG:                    {png_path}")

def main() -> None:
    code = BCH(TOTAL_BITS, BIT_CORRECTIONS)
    print(code.summary())

    if RUN_SELF_TEST:
        passed, total = self_test(code, SELF_TEST_TRIALS, RANDOM_SEED)
        print(f"\nSelf-test: {passed}/{total} passed")
        if passed != total:
            raise SystemExit("Self-test failed; refusing to process data.")

    if not LOAD_FILE_PATHS:
        print(
            "\nLOAD_FILE_PATHS is empty, so no data files were processed.\n"
            "Add one or more paths at the top of this script and run it again."
        )
        return

    for file_path in LOAD_FILE_PATHS:
        process_file(file_path, code)


if __name__ == "__main__":
    main()
