#!/usr/bin/env python3
import re
import sys
from pathlib import Path

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

# strip ANSI escapes
ansi_escape = re.compile(r'\x1b\[[0-9;]*[mK]')
def clean_ansi(txt: str) -> str:
    return ansi_escape.sub('', txt)

def float_to_mathematica(v: float) -> str:
    """
    Convert a Python float to Mathematica scientific or decimal notation.
    E.g. 4.1e-06 -> "4.1*^-06", -0.02 -> "-0.02".
    """
    s = f"{v:.6e}"
    m_str, e_str = s.split('e')
    e_val = int(e_str)
    m_str = m_str.rstrip('0').rstrip('.')
    if e_val == 0:
        return repr(v)
    sign = '+' if e_val > 0 else '-'
    exp_str = f"{abs(e_val):02d}"
    return f"{m_str}*^{sign}{exp_str}"

# -----------------------------------------------------------------------------
# Main parsing
# -----------------------------------------------------------------------------

def parse_table(filename: str):
    # regex for "number +- error"
    num_re = re.compile(
        r"\s*([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)"  # central
        r"\s*\+\-\s*"
        r"([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)"      # error
    )

    text = Path(filename).read_text()
    lines = clean_ansi(text).splitlines()

    # find the two '+'–separator lines around the header
    sep_idxs = [i for i, L in enumerate(lines) if L.startswith('+') and L.endswith('+')]
    if len(sep_idxs) < 2:
        raise RuntimeError("Couldn't find table separators")
    header_sep = sep_idxs[1]
    header_line = lines[header_sep - 1]
    sep_line    = lines[header_sep]

    # compute column‐slice boundaries from '+' chars
    boundaries = [i for i, ch in enumerate(sep_line) if ch == '+']
    # extract header names
    headers = [
        header_line[start+1:end].strip()
        for start, end in zip(boundaries, boundaries[1:])
    ]

    data = {}
    for L in lines[header_sep+1:]:
        if L.startswith('+'):
            break
        if not L.startswith('|'):
            continue

        # slice the row into cells
        cells = [
            L[start+1:end].strip()
            for start, end in zip(boundaries, boundaries[1:])
        ]
        row_id = cells[0]
        # skip only the summary "Sum" row, but keep "Total"
        if row_id == "Sum":
            continue

        row_map = {}

        # --- Column 1: measurement (+ inline metrics) ---
        raw = cells[1]
        if '|' in raw:
            meas_part, metrics_part = raw.split('|', 1)
        else:
            meas_part, metrics_part = raw, ''
        m = num_re.match(meas_part)
        if not m:
            raise ValueError(f"Can't parse measurement in: {meas_part!r}")
        central = float(m.group(1))
        error   = float(m.group(2))

        # derive the Mathematica key from the header
        key_full = re.search(r"\((.*?)\)", headers[1]).group(1)
        suffix   = key_full.replace('mu_r_x_', '').replace('_', '')
        row_map[f"mux{suffix}"] = (central, error)

        # pull out any k=v metrics in that same cell
        for kv in metrics_part.split('|'):
            if '=' not in kv:
                continue
            k, v = kv.split('=', 1)
            v = v.strip()
            if re.fullmatch(r'[+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?', v):
                row_map[k] = float(v)
            else:
                row_map[k] = v

        # --- Any additional “plain” measurement columns ---
        # (this will be empty if you only ever have two columns)
        for idx in range(2, len(cells)):
            cell = cells[idx]
            m2 = num_re.match(cell)
            if not m2:
                raise ValueError(f"Can't parse {headers[idx] if idx < len(headers) else 'col'+str(idx)}: {cell!r}")
            cval = float(m2.group(1))
            err  = float(m2.group(2))

            if idx < len(headers):
                key_full = re.search(r"\((.*?)\)", headers[idx]).group(1)
                suffix   = key_full.replace('mu_r_x_', '').replace('_','')
            else:
                suffix = str(idx)
            row_map[f"mux{suffix}"] = (cval, err)

        data[row_id] = row_map

    return data

def format_mathematica(data):
    rows = []
    for row_id, fmap in data.items():
        parts = []
        for k, v in fmap.items():
            if isinstance(v, tuple):
                c_str = float_to_mathematica(v[0])
                e_str = float_to_mathematica(v[1])
                parts.append(f'"{k}" -> {{{c_str}, {e_str}}}')
            elif isinstance(v, float):
                parts.append(f'"{k}" -> {float_to_mathematica(v)}')
            else:
                parts.append(f'"{k}" -> "{v}"')
        rows.append(f'"{row_id}" -> <|{", ".join(parts)}|>')
    return '<|' + ', '.join(rows) + '|>'

# -----------------------------------------------------------------------------
# CLI entrypoint
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <table-file>", file=sys.stderr)
        sys.exit(1)
    data = parse_table(sys.argv[1])
    print(format_mathematica(data))

