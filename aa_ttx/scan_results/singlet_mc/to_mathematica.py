#!/usr/bin/env python3
import re
import sys
import glob

import re, pathlib

ansi_escape = re.compile(r'\x1b\[[0-9;]*[mK]')
def clean_ansi(txt):
    return ansi_escape.sub('', txt)

def float_to_mathematica(v):
    """
    Convert a Python float to Mathematica scientific or decimal notation.
    E.g. 4.1e-06 -> "4.1*^-06", -0.02 -> "-0.02".
    """
    # Use scientific format with 6-digit mantissa
    s = f"{v:.6e}"
    m_str, e_str = s.split('e')
    e_val = int(e_str)
    # Strip trailing zeros and dot from mantissa
    m_str = m_str.rstrip('0').rstrip('.')
    if e_val == 0:
        # Return plain Python repr for exact decimal
        return repr(v)
    # Prepare exponent with sign and two digits
    sign = '+' if e_val > 0 else '-'
    exp = abs(e_val)
    exp_str = f"{exp:02d}"
    return f"{m_str}*^{sign}{exp_str}"


def parse_table(filename):
    num_re = re.compile(
        r"\s*([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)"  # central value
        r"\s*\+\-\s*"
        r"([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)"      # error
    )

    with open(filename) as f:
        lines = clean_ansi(f.read()).splitlines()

    # Locate header separators (lines starting and ending with '+')
    sep_idxs = [i for i, l in enumerate(lines) if l.startswith('+') and l.endswith('+')]
    if len(sep_idxs) < 2:
        raise RuntimeError("Couldn't find table separators")
    header_sep = sep_idxs[1]
    header_line = lines[header_sep - 1]
    sep_line    = lines[header_sep]

    # Determine column boundaries from '+' positions in the separator line
    boundaries = [i for i, ch in enumerate(sep_line) if ch == '+']
    headers = [
        header_line[start+1:end].strip()
        for start, end in zip(boundaries, boundaries[1:])
    ]

    data = {}
    for line in lines[header_sep + 1:]:
        if line.startswith('+'):
            break
        if not line.startswith('|'):
            continue

        # Slice out each cell by the same boundaries
        cells = [
            line[start+1:end].strip()
            for start, end in zip(boundaries, boundaries[1:])
        ]
        row_id = cells[0]
        if row_id == "Sum":
            continue

        row_map = {}

        # Column 1: measurement + metrics
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

        # Build the key name from the header, e.g. "(mu_r_x_1)" → "mux1"
        key_full = re.search(r"\((.*?)\)", headers[1]).group(1)
        suffix   = key_full.replace('mu_r_x_', '').replace('_', '')
        key1     = f"mux{suffix}"
        row_map[key1] = (central, error)

        # Metrics after the first '|'
        for kv in metrics_part.split('|'):
            if '=' not in kv:
                continue
            k, v = kv.split('=', 1)
            v = v.strip()
            # convert numeric strings to float
            if re.fullmatch(r'[+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?', v):
                row_map[k] = float(v)
            else:
                row_map[k] = v

        # Columns 2 & 3: plain measurements (only if present)
        for idx in (2, 3):
            if idx >= len(headers) or idx >= len(cells):
                continue
            cell = cells[idx]
            if not cell:
                continue
            m2 = num_re.match(cell)
            if not m2:
                raise ValueError(f"Can't parse {headers[idx]}: {cell!r}")
            cval = float(m2.group(1))
            err  = float(m2.group(2))
            key_full = re.search(r"\((.*?)\)", headers[idx]).group(1)
            suffix   = key_full.replace('mu_r_x_', '').replace('_', '')
            key_name = f"mux{suffix}"
            row_map[key_name] = (cval, err)

        data[row_id] = row_map

    return data

def parse_table_OLD(filename):
    num_re = re.compile(
        r"\s*([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)"  # central
        r"\s*\+\-\s*"
        r"([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?)"      # error
    )

    with open(filename) as f:
        lines = clean_ansi(f.read()).splitlines()

    # Locate header separator
    sep_idxs = [i for i,l in enumerate(lines) if l.startswith('+') and l.endswith('+')]
    if len(sep_idxs) < 2:
        raise RuntimeError("Couldn't find table separators")
    header_sep = sep_idxs[1]
    header_line = lines[header_sep-1]
    sep_line    = lines[header_sep]

    # Determine column boundaries from '+' positions
    boundaries = [i for i,ch in enumerate(sep_line) if ch == '+']
    headers = []
    for start, end in zip(boundaries, boundaries[1:]):
        headers.append(header_line[start+1:end].strip())

    data = {}
    for line in lines[header_sep+1:]:
        if line.startswith('+'):
            break
        if not line.startswith('|'):
            continue

        cells = [line[start+1:end].strip() for start, end in zip(boundaries, boundaries[1:])]
        row_id = cells[0]
        if row_id == "Sum":
            continue

        row_map = {}
        # Column 1: measurement + metrics
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
        key_full = re.search(r"\((.*?)\)", headers[1]).group(1)
        suffix   = key_full.replace('mu_r_x_', '').replace('_','')
        key1     = f"mux{suffix}"
        row_map[key1] = (central, error)

        # Metrics
        for kv in metrics_part.split('|'):
            if '=' not in kv:
                continue
            k, v = kv.split('=', 1)
            v = v.strip()
            # convert numeric strings to float, else keep
            if re.fullmatch(r'[+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?', v):
                row_map[k] = float(v)
            else:
                row_map[k] = v

        # Columns 2 & 3: plain measurements
        for idx in (2, 3):
            print(cells, idx)
            cell = cells[idx]
            m2 = num_re.match(cell)
            if not m2:
                raise ValueError(f"Can't parse {headers[idx]}: {cell!r}")
            cval = float(m2.group(1))
            err  = float(m2.group(2))
            key_full = re.search(r"\((.*?)\)", headers[idx]).group(1)
            suffix   = key_full.replace('mu_r_x_', '').replace('_','')
            key_name = f"mux{suffix}"
            row_map[key_name] = (cval, err)

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

id_extractor = re.compile(r"^nnlo_all_singlets_(\d+)(_?.*)_sqrts_(\d+)_(\d+)_(.*)\.txt$")
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <table-file>", file=sys.stderr)
        sys.exit(1)
    all_file_paths = glob.glob(sys.argv[1])
    all_files = []
    for fname in  glob.glob(sys.argv[1]):
        run_id = int(id_extractor.match(fname).group(1))
        sqrt_s = float(id_extractor.match(fname).group(3)+"."+id_extractor.match(fname).group(4))
        all_files.append((run_id, sqrt_s, fname))
    all_files = sorted(all_files, key=lambda k: k[1])
    all_res = []
    for i_scan, sqrt_s, f in all_files:
        print("Processing file '%s'..."%f)
        data = parse_table(f)
        all_res.append("{ %f, %s }"%(sqrt_s, format_mathematica(data)))
    with open('scan_results_mm.txt','w') as f:
        f.write("{\n"+",\n".join(all_res)+"\n}")
