#!/usr/bin/env python3
"""
21_yn00Harvester_streamed_v5_3.py

Streamed harvester for PAML yn00 outputs — v5.3 (non-invasive codon-pos normalization).

Non-invasive edits from v5/v5.2:
 - codon_pos_base_seq: per-sequence entries no longer include an "avg" slot.
   Each per-sequence entry has seq_index and seq_name set to the same label (e.g. "seq054").
   If the file uses "# 1: NAME" lines, those are normalized to "seq001", "seq002", ...
 - codon_pos_base_avg: stores the Average block as a JSON *list containing a single dict*:
     [{ "seq_index":"Average", "seq_name":"Average", "pos1":[...], "pos2":[...], "pos3":[...] }]
 - No other parsing/column-writing logic changed.
"""
from __future__ import annotations
import os, re, glob, csv, json
from pathlib import Path
from typing import Tuple, List, Dict, Any

# tolerant numeric regex (floats, sci, inf, nan)
NUM_RE = re.compile(
    r'(?:[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?|inf|infty|-inf|nan|-nan)',
    re.IGNORECASE
)

# -------------------------
# Existing fields (unchanged) followed by added/new fields appended at the end
# -------------------------
FIELDNAMES = [
    # existing core fields (unchanged order & names)
    "gene_id","seq1","seq2",
    "S","N","t","kappa","omega","dN","dN_SE","dS","dS_SE",
    "L_i","Ns_i","Nv_i","A_i","B_i",
    "LWL85_dS","LWL85_dN","LWL85_omega","LWL85_S","LWL85_N",
    "LWL85m_dS","LWL85m_dN","LWL85m_omega","LWL85m_S","LWL85m_N","LWL85m_rho",
    "LPB93_dS","LPB93_dN","LPB93_omega",
    "seq_count",
    "PreGapFilter_nSites","PostGapFilter_nSites","ns","ls",
    # ---- newly requested columns (appended) ----
    "matrix_dS","matrix_dN","matrix_t",
    "seq_labels",
    "codon_counts",
    "GC_pos1","GC_pos2","GC_pos3","GC_total",
    "pair_n_identical_codons","pair_n_gapped_codons","pair_n_ambiguous_codons","pair_pct_identity",
    "site_pattern_counts",
    "dS_is_zero","any_nan_inf","omega_placeholder_flag",
    # NEW (from FreeRatio logic) - non-invasive additions
    "codon_usage_counts","codon_pos_base_avg","codon_pos_base_seq"
]

# backward mapping used by writer to preserve internal keys
BACKWARD_SCALAR_MAP = {
    "S": "YN00_S", "N": "YN00_N", "t": "YN00_t", "kappa": "YN00_kappa",
    "omega": "YN00_omega", "dN": "YN00_dN", "dN_SE": "YN00_dN_SE",
    "dS": "YN00_dS", "dS_SE": "YN00_dS_SE",
    "LWL85_dS": "LWL85_dS", "LWL85_dN": "LWL85_dN", "LWL85_omega": "LWL85_omega",
    "LWL85_S": "LWL85_S", "LWL85_N": "LWL85_N",
    "LWL85m_dS": "LWL85m_dS", "LWL85m_dN": "LWL85m_dN", "LWL85m_omega": "LWL85m_omega",
    "LWL85m_S": "LWL85m_S", "LWL85m_N": "LWL85m_N", "LWL85m_rho": "LWL85m_rho",
    "LPB93_dS": "LPB93_dS", "LPB93_dN": "LPB93_dN", "LPB93_omega": "LPB93_omega",
    "seq_count": "seq_count",
    # Pre/Post/ ns/ls adopted into parser rows as keys; writer reads them directly
}

# -------------------------
# Additional regexes and helpers (from FreeRatio harvester) for codon blocks
# -------------------------
SEQNUM_RE        = re.compile(r'^\s*#\s*([0-9]+)\s*:\s*(.*)$')
CODON_POS_HEADER_RE = re.compile(r'Codon position x base', re.IGNORECASE)
POS_BASE_RE         = re.compile(r'position\s*(\d+)\s*:\s*(.*)', re.IGNORECASE)
CODON_SUMS_HEADER_RE = re.compile(r'Sums of codon usage counts', re.IGNORECASE)
CODON_COUNT_RE      = re.compile(r'\b([ACGT]{3})\b\s*([0-9]+)', re.IGNORECASE)

# -------------------------
# Helper functions
# -------------------------
def safe_float_token(tok: str):
    if tok is None:
        return None
    t = str(tok).strip().lower()
    if t in ("inf","infty","+inf","-inf","nan","-nan"):
        return None
    try:
        return float(tok)
    except:
        return None

def has_zero_sites_after_gaps(text: str) -> bool:
    return bool(re.search(r'after\s+deleting\s+gaps[^\d\n]*\b0\b\s+sites', text, re.IGNORECASE))

def build_seq_token_list(text: str) -> List[str]:
    tokens: List[str] = []
    for m in re.finditer(r'^\s*(seq\d+)\b', text, flags=re.I | re.M):
        lbl = m.group(1).strip()
        if lbl not in tokens:
            tokens.append(lbl)
    return tokens

def map_to_seq_label(token: str, seq_tokens: List[str]) -> str:
    tok = token.strip()
    if re.fullmatch(r'seq\d+', tok, re.I):
        return tok
    if re.fullmatch(r'\d+', tok):
        idx = int(tok)
        if 1 <= idx <= len(seq_tokens):
            return seq_tokens[idx-1]
        return f"seq{idx:03d}"
    return tok

# -------------------------
# New extractors (non-invasive; improved codon-pos parsing for yn00 outputs)
# -------------------------
def extract_codon_usage_counts_dict(text: str) -> str:
    m = re.search(r'\bSums\b(.*?)(?:\n\n|\Z)', text, re.S | re.I)
    block = m.group(1) if m else ""
    if not block:
        m2 = re.search(r'Codon usage for each species(.*?)(?:\n\n|\Z)', text, re.S | re.I)
        block = m2.group(1) if m2 else ""
    if not block:
        return ""
    counts: Dict[str,int] = {}
    for codon, cnt in CODON_COUNT_RE.findall(block):
        try:
            counts[codon.upper()] = counts.get(codon.upper(), 0) + int(cnt)
        except:
            pass
    if not counts:
        return ""
    return json.dumps(counts, separators=(',', ':'))

def extract_codon_pos_base_seq_and_avg(text: str) -> Tuple[str, str]:
    """
    Robust extractor for yn00 codon-position block.

    Returns (codon_pos_list_json, avg_json) where:
      - codon_pos_list_json is JSON list of per-sequence dicts:
           { 'seq_index': 'seq054', 'seq_name': 'seq054',
             'pos1':[T,C,A,G], 'pos2':[...], 'pos3':[...] }
      - avg_json is JSON list with a single dict:
           [{ 'seq_index': 'Average', 'seq_name': 'Average', 'pos1':[...], 'pos2':[...], 'pos3':[...] }]
      - If a per-seq block lacks any pos values it will be omitted or have nulls for those positions.
    """
    lines = text.splitlines()
    start_idx = None
    for i, ln in enumerate(lines):
        if CODON_POS_HEADER_RE.search(ln):
            start_idx = i
            break
    if start_idx is None:
        return "", ""

    seq_entries: List[Dict[str, Any]] = []
    avg_block_values: Dict[str, List[float]] = {}
    curr: Dict[str, Any] | None = None

    i = start_idx + 1
    while i < len(lines):
        raw = lines[i]
        line = raw.rstrip("\n")
        stripped = line.strip()

        # sequence header forms:
        m_seqnum = SEQNUM_RE.match(line)            # "# 1: name"
        m_seqlabel = re.match(r'^\s*(seq\d+)\b\s*(.*)$', line, re.I)  # "seq052    " or "seq052    rest"

        if m_seqnum:
            # commit previous
            if curr is not None:
                seq_entries.append(curr)
            try:
                idx_val = int(m_seqnum.group(1))
            except:
                idx_val = None
            # normalize label to seq### so seq_index and seq_name are identical
            label = f"seq{int(m_seqnum.group(1)):03d}" if idx_val is not None else m_seqnum.group(1).strip()
            curr = {"seq_index": label, "seq_name": label, "pos1": None, "pos2": None, "pos3": None}
            i += 1
            continue

        elif m_seqlabel and not stripped.lower().startswith("position"):
            # commit previous
            if curr is not None:
                seq_entries.append(curr)
            label = m_seqlabel.group(1).strip()
            # set both seq_index and seq_name to the same label
            curr = {"seq_index": label, "seq_name": label, "pos1": None, "pos2": None, "pos3": None}
            i += 1
            continue

        # Average block start
        if re.match(r'^\s*Average\b', line, re.I):
            # parse subsequent position lines for avg
            j = i + 1
            avg_temp: Dict[int, List[float]] = {}
            while j < len(lines):
                lm = POS_BASE_RE.match(lines[j])
                if not lm:
                    break
                posnum = int(lm.group(1))
                rest = lm.group(2)
                base_map = {}
                for bm in re.finditer(r'(T|C|A|G)\s*:\s*([-\d\.Ee+]+)', rest, re.I):
                    base = bm.group(1).upper()
                    val = safe_float_token(bm.group(2))
                    base_map[base] = val
                ordered = [ base_map.get(b) for b in ("T","C","A","G") ]
                avg_temp[posnum] = ordered
                j += 1
            # populate avg_block_values
            if avg_temp:
                for p in (1,2,3):
                    vals = avg_temp.get(p)
                    if vals:
                        avg_block_values[f"pos{p}"] = vals
            i = j
            # commit any open curr (we treat Average separately)
            if curr is not None:
                seq_entries.append(curr)
                curr = None
            continue

        # position lines for current sequence
        pm = POS_BASE_RE.match(line)
        if pm and curr is not None:
            posnum = int(pm.group(1))
            rest = pm.group(2)
            base_map = {}
            for bm in re.finditer(r'(T|C|A|G)\s*:\s*([-\d\.Ee+]+)', rest, re.I):
                base = bm.group(1).upper()
                val = safe_float_token(bm.group(2))
                base_map[base] = val
            ordered = [ base_map.get(b) for b in ("T","C","A","G") ]
            if posnum == 1:
                curr["pos1"] = ordered
            elif posnum == 2:
                curr["pos2"] = ordered
            elif posnum == 3:
                curr["pos3"] = ordered
            i += 1
            continue

        # blank line: commit current if present
        if not stripped:
            if curr is not None:
                seq_entries.append(curr)
                curr = None
            i += 1
            continue

        # fallback: move on
        i += 1

    # final commit
    if curr is not None:
        seq_entries.append(curr)

    # Build JSON outputs
    # Filter out totally empty entries (optional)
    filtered_seq_entries = []
    for e in seq_entries:
        # keep at least if any pos present
        if any(e.get(f"pos{p}") is not None for p in (1,2,3)):
            filtered_seq_entries.append(e)
    codon_pos_json = json.dumps(filtered_seq_entries, separators=(',', ':')) if filtered_seq_entries else ""

    # avg JSON: wrap into a list with a labeled Average entry if present
    avg_json = ""
    if avg_block_values:
        avg_entry = {
            "seq_index": "Average",
            "seq_name": "Average",
            "pos1": avg_block_values.get("pos1"),
            "pos2": avg_block_values.get("pos2"),
            "pos3": avg_block_values.get("pos3")
        }
        avg_json = json.dumps([avg_entry], separators=(',', ':'))

    return codon_pos_json, avg_json

# -------------------------
# Reused extractors from v5 for other fields (unchanged)
# -------------------------
def extract_pregap_nsites(text: str) -> int | None:
    m = re.search(r'Before\s+deleting\s+alignment\s+gaps\s*\n\s*([0-9]+)\s+([0-9]+)', text, re.I)
    if m:
        try:
            return int(m.group(2))
        except:
            return None
    return None

def extract_postgap_nsites(text: str) -> int | None:
    m = re.search(r'After\s+deleting\s+gaps\.\s*([0-9]+)\s+sites', text, re.I)
    if m:
        try:
            return int(m.group(1))
        except:
            return None
    return None

def extract_ns_ls(text: str) -> Tuple[int|None, int|None]:
    m = re.search(r'ns\s*=\s*([0-9]+)\s+ls\s*=\s*([0-9]+)', text, re.I)
    if m:
        try:
            return int(m.group(1)), int(m.group(2))
        except:
            return None, None
    return None, None

def extract_seq_labels_json(text: str) -> str:
    seqs = build_seq_token_list(text)
    if not seqs:
        return ""
    return json.dumps(seqs)

def extract_site_pattern_counts(text: str) -> str:
    m = re.search(r'Printing out site pattern counts\s*(.*?)\n\n', text, re.S | re.I)
    if not m:
        m2 = re.search(r'Printing out site pattern counts\s*(.*?)(?:\n[A-Za-z]|\Z)', text, re.S | re.I)
        if not m2:
            return ""
        sblock = m2.group(1)
    else:
        sblock = m.group(1)
    nums = re.findall(r'\b\d+\b', sblock)
    if not nums:
        return ""
    return "[" + ", ".join(nums) + "]"

def extract_codon_counts_json(text: str) -> str:
    m = re.search(r'\bSums\b(.*?)(?:\n\n|\Z)', text, re.S | re.I)
    block = m.group(1) if m else ""
    if not block:
        m2 = re.search(r'Codon usage for each species(.*?)(?:\n\n|\Z)', text, re.S | re.I)
        block = m2.group(1) if m2 else ""
    if not block:
        return ""
    counts = {}
    for codon, cnt in re.findall(r'\b([ACGT]{3})\b\s+([0-9]+)', block, re.I):
        counts[codon.upper()] = counts.get(codon.upper(), 0) + int(cnt)
    if not counts:
        return ""
    return json.dumps(counts, separators=(',', ':'))

def extract_codon_pos_gc(text: str) -> Tuple[float|None, float|None, float|None]:
    m = re.search(r'Average\s*(.*?)\n\n', text, re.S | re.I)
    snippet = m.group(1) if m else text
    pos_vals = {}
    for pos in (1,2,3):
        pm = re.search(r'position\s+' + str(pos) + r'\s*:\s*(.*)', text, re.I)
        if pm:
            vals = {}
            for base_match in re.finditer(r'(T|C|A|G)\s*:\s*([0-9]*\.?[0-9]+)', pm.group(1), re.I):
                base = base_match.group(1).upper()
                val = float(base_match.group(2))
                vals[base] = val
            pos_vals[pos] = vals.get('G') + vals.get('C') if ('G' in vals and 'C' in vals) else None
        else:
            pos_vals[pos] = None
    g1 = pos_vals.get(1)
    g2 = pos_vals.get(2)
    g3 = pos_vals.get(3)
    if any(v is not None for v in (g1,g2,g3)):
        nums = [v for v in (g1,g2,g3) if v is not None]
        gc_total = sum(nums)/len(nums) if nums else None
        return g1, g2, g3, gc_total
    return None, None, None, None

# -----------------------------------------
# Robust extractor for "Before deleting alignment gaps" block (copied)
# -----------------------------------------
def extract_before_alignment_seqs(text: str) -> Dict[str, List[str]]:
    seqmap: Dict[str, List[str]] = {}

    start_m = re.search(r'Before\s+deleting\s+alignment\s+gaps', text, re.I)
    if not start_m:
        return {}

    start_pos = start_m.end()
    after_m = re.search(r'After\s+deleting\s+gaps', text[start_pos:], re.I)
    end_pos = start_pos + after_m.start() if after_m else None

    if end_pos:
        block = text[start_pos:end_pos]
    else:
        lines = text[start_pos:].splitlines()
        capture_lines = []
        for ln in lines:
            capture_lines.append(ln)
            if len(capture_lines) > 2000:
                break
        block = "\n".join(capture_lines)

    last_label = None
    for line in block.splitlines():
        if not line.strip():
            last_label = None
            continue

        m = re.match(r'^\s*(seq\d+)\b\s*(.*)$', line, re.I)
        if m:
            label = m.group(1).strip()
            rest = m.group(2).strip()
            toks = [t for t in re.split(r'\s+', rest) if t != '']
            seqmap.setdefault(label, []).extend(toks)
            last_label = label
        else:
            if last_label is not None:
                rest = line.strip()
                toks = [t for t in re.split(r'\s+', rest) if t != '']
                seqmap.setdefault(last_label, []).extend(toks)

    for k, v in list(seqmap.items()):
        seqmap[k] = [t.upper() for t in v if t != '']

    return seqmap

# -----------------------------------------
# Robust pairwise codon stats computation (copied)
# -----------------------------------------
def compute_pair_codons_stats(seqmap: Dict[str, List[str]], s1: str, s2: str) -> Tuple[int,int,int,float]:
    toks1 = seqmap.get(s1, []) or seqmap.get(s1.upper(), [])
    toks2 = seqmap.get(s2, []) or seqmap.get(s2.upper(), [])
    if not toks1 or not toks2:
        return 0, 0, 0, None

    n = min(len(toks1), len(toks2))
    identical = 0
    gapped = 0
    ambiguous = 0

    def is_gap(tok: str) -> bool:
        return ('-' in tok) or (tok == '---')

    def is_valid_codon(tok: str) -> bool:
        return bool(re.fullmatch(r'[ACGT]{3}', tok, re.I))

    def is_ambiguous(tok: str) -> bool:
        if tok == '':
            return True
        if is_gap(tok):
            return False
        return not is_valid_codon(tok)

    for i in range(n):
        a = toks1[i].upper()
        b = toks2[i].upper()
        if is_gap(a) or is_gap(b):
            gapped += 1
            continue
        if is_ambiguous(a) or is_ambiguous(b):
            ambiguous += 1
            continue
        if a == b:
            identical += 1

    denom = n - gapped - ambiguous
    pct = (identical / denom * 100.0) if denom and denom > 0 else None
    return identical, gapped, ambiguous, pct

def parse_matrix_file_into_lower_tri_values(matrix_path: Path) -> List[float]:
    try:
        txt = matrix_path.read_text(encoding='utf-8', errors='replace')
    except:
        return []
    nums = NUM_RE.findall(txt)
    vals = []
    for n in nums:
        try:
            vals.append(float(n))
        except:
            pass
    return vals

# -------------------------
# Main parser (keeps prior behavior) - returns rows and message
# -------------------------
def parse_yn00_file_rows(filepath: str) -> Tuple[List[Dict[str,Any]], str]:
    path = Path(filepath)
    try:
        gene_name = Path(filepath).parts[-3]
    except Exception:
        gene_name = path.parent.parent.name

    try:
        content = path.read_text(encoding='utf-8', errors='replace')
    except Exception as e:
        return [], f"⚠ {gene_name}: Read Error - {e}"

    # gap check
    if has_zero_sites_after_gaps(content):
        return [], f"⚠ {gene_name}: Gap Warning-- 0 sites after gap removal."

    # sanity
    if "Time used:" not in content and "Yang & Nielsen" not in content:
        return [], f"⚠ {gene_name}: Parsed 0 branches (File incomplete)"

    # extract repeating values early
    pre_n = extract_pregap_nsites(content)
    post_n = extract_postgap_nsites(content)
    ns_val, ls_val = extract_ns_ls(content)
    seq_tokens = build_seq_token_list(content)
    seq_count = len(seq_tokens)
    seq_labels_json = extract_seq_labels_json(content)
    site_pattern_vec = extract_site_pattern_counts(content)
    codon_counts_json = extract_codon_counts_json(content)

    # codon pos GC extraction (returns 4-tuple)
    g1, g2, g3, gtotal = extract_codon_pos_gc(content)

    # pre-alignment seq map for pairwise counts (before deletion)
    before_seqmap = extract_before_alignment_seqs(content)

    # ALSO: extract the FreeRatio-like codon blocks into JSON strings (non-invasive)
    codon_usage_counts_json = extract_codon_usage_counts_dict(content)  # aggregated codon counts
    codon_pos_base_seq_json, codon_pos_base_avg_json = extract_codon_pos_base_seq_and_avg(content)  # per-seq codon pos base table + average

    # merged pairwise rows (keep original behavior)
    merged: Dict[tuple, Dict[str,Any]] = {}

    yn_master_re = re.compile(
        r"^\s*(\S+)\s+(\S+)\s+"
        r"(" + NUM_RE.pattern + r")\s+(" + NUM_RE.pattern + r")\s+"
        r"(" + NUM_RE.pattern + r")\s+(" + NUM_RE.pattern + r")\s+"
        r"(" + NUM_RE.pattern + r")\s+"
        r"(" + NUM_RE.pattern + r")\s+\+-\s+(" + NUM_RE.pattern + r")\s+"
        r"(" + NUM_RE.pattern + r")\s+\+-\s+(" + NUM_RE.pattern + r")",
        re.I | re.M
    )

    for m in yn_master_re.finditer(content):
        raw_s1, raw_s2 = m.group(1), m.group(2)
        s1 = map_to_seq_label(raw_s1, seq_tokens)
        s2 = map_to_seq_label(raw_s2, seq_tokens)
        key = tuple(sorted([s1, s2]))
        merged[key] = {
            "gene_id": gene_name, "seq1": s1, "seq2": s2,
            # keep internal keys as before
            "YN00_S": safe_float_token(m.group(3)),
            "YN00_N": safe_float_token(m.group(4)),
            "YN00_t": safe_float_token(m.group(5)),
            "YN00_kappa": safe_float_token(m.group(6)),
            "YN00_omega": safe_float_token(m.group(7)),
            "YN00_dN": safe_float_token(m.group(8)),
            "YN00_dN_SE": safe_float_token(m.group(9)),
            "YN00_dS": safe_float_token(m.group(10)),
            "YN00_dS_SE": safe_float_token(m.group(11)),
            "seq_count": seq_count,
            # repeating metadata
            "PreGapFilter_nSites": pre_n,
            "PostGapFilter_nSites": post_n,
            "ns": ns_val,
            "ls": ls_val,
            # extra enriched fields (fill below)
            "seq_labels": seq_labels_json,
            "site_pattern_counts": site_pattern_vec,
            "codon_counts": codon_counts_json,
            "GC_pos1": g1, "GC_pos2": g2, "GC_pos3": g3, "GC_total": gtotal,
            # NEW: insert codon summaries collected above
            "codon_usage_counts": codon_usage_counts_json,
            "codon_pos_base_seq": codon_pos_base_seq_json,
            "codon_pos_base_avg": codon_pos_base_avg_json
        }

    # detailed blocks (unchanged behavior)
    block_re = re.compile(
        r"^\s*(\d+)\s*\(\s*([^)]+?)\s*\)\s+vs\.\s+(\d+)\s*\(\s*([^)]+?)\s*\)(.*?)(?=\n\s*\d+\s*\(|\n\s*\(C\)\s+LWL85|\Z)",
        re.I | re.M | re.S
    )

    def parse_triple_sum(label: str, block_text: str):
        m = re.search(rf"{re.escape(label)}\(i\)\s*:\s*(.*)", block_text, re.I)
        if not m:
            return None
        tail = m.group(1)
        nums = NUM_RE.findall(tail)
        out: Dict[str, Any] = {}
        if len(nums) >= 1:
            out[f"{label}_1"] = safe_float_token(nums[0]) if len(nums) >= 1 else None
            out[f"{label}_2"] = safe_float_token(nums[1]) if len(nums) >= 2 else None
            out[f"{label}_3"] = safe_float_token(nums[2]) if len(nums) >= 3 else None
        sum_m = re.search(r"sum\s*=\s*(" + NUM_RE.pattern + r")", tail, re.I)
        if sum_m:
            out[f"{label}_sum"] = safe_float_token(sum_m.group(1))
        return out or None

    for bm in block_re.finditer(content):
        idx1, name1 = bm.group(1).strip(), bm.group(2).strip()
        idx2, name2 = bm.group(3).strip(), bm.group(4).strip()
        block_text = bm.group(5)
        s1 = name1 if name1.lower().startswith("seq") else map_to_seq_label(name1, seq_tokens)
        s2 = name2 if name2.lower().startswith("seq") else map_to_seq_label(name2, seq_tokens)

        candidates = [
            tuple(sorted([s1, s2])),
            tuple(sorted([map_to_seq_label(idx1, seq_tokens), map_to_seq_label(idx2, seq_tokens)])),
            tuple(sorted([name1, name2])),
            tuple(sorted([f"seq{int(idx1):03d}", f"seq{int(idx2):03d}"])) if idx1.isdigit() and idx2.isdigit() else None
        ]
        mapped_key = None
        for c in candidates:
            if not c:
                continue
            if c in merged:
                mapped_key = c
                break

        if mapped_key is None:
            s1_label = map_to_seq_label(idx1, seq_tokens) if idx1.isdigit() else s1
            s2_label = map_to_seq_label(idx2, seq_tokens) if idx2.isdigit() else s2
            mapped_key = tuple(sorted([s1_label, s2_label]))
            if mapped_key not in merged:
                merged[mapped_key] = {
                    "gene_id": gene_name, "seq1": s1_label, "seq2": s2_label,
                    "YN00_S": None, "YN00_N": None, "YN00_t": None, "YN00_kappa": None,
                    "YN00_omega": None, "YN00_dN": None, "YN00_dN_SE": None,
                    "YN00_dS": None, "YN00_dS_SE": None,
                    "seq_count": seq_count,
                    "PreGapFilter_nSites": pre_n,
                    "PostGapFilter_nSites": post_n,
                    "ns": ns_val,
                    "ls": ls_val,
                    "seq_labels": seq_labels_json,
                    "site_pattern_counts": site_pattern_vec,
                    "codon_counts": codon_counts_json,
                    "GC_pos1": g1, "GC_pos2": g2, "GC_pos3": g3, "GC_total": gtotal,
                    "codon_usage_counts": codon_usage_counts_json,
                    "codon_pos_base_seq": codon_pos_base_seq_json,
                    "codon_pos_base_avg": codon_pos_base_avg_json
                }

        row = merged[mapped_key]

        for label in ("L", "Ns", "Nv"):
            parsed = parse_triple_sum(label, block_text)
            if parsed:
                row.update(parsed)

        # A(i)/B(i)
        for label in ("A", "B"):
            m = re.search(rf"{re.escape(label)}\(i\)\s*:\s*(.*)", block_text, re.I)
            if m:
                nums = NUM_RE.findall(m.group(1))
                for idx in range(min(3, len(nums))):
                    row[f"{label}_{idx+1}"] = safe_float_token(nums[idx])

        # LWL85 / LWL85m / LPB93 parsing (unchanged behavior)
        lwl85 = re.search(
            r"LWL85\s*:\s*dS\s*=\s*(" + NUM_RE.pattern + r")\s*dN\s*=\s*(" + NUM_RE.pattern +
            r")\s*w\s*=\s*(" + NUM_RE.pattern + r")\s*S\s*=\s*(" + NUM_RE.pattern + r")\s*N\s*=\s*(" + NUM_RE.pattern + r")",
            block_text, re.I
        )
        if lwl85:
            vals = lwl85.groups()
            row["LWL85_dS"] = safe_float_token(vals[0])
            row["LWL85_dN"] = safe_float_token(vals[1])
            row["LWL85_omega"] = safe_float_token(vals[2])
            row["LWL85_S"] = safe_float_token(vals[3])
            row["LWL85_N"] = safe_float_token(vals[4])

        lwl85m = re.search(
            r"LWL85m\s*:\s*dS\s*=\s*(" + NUM_RE.pattern + r")\s*dN\s*=\s*(" + NUM_RE.pattern +
            r")\s*w\s*=\s*(" + NUM_RE.pattern + r")\s*S\s*=\s*(" + NUM_RE.pattern + r")\s*N\s*=\s*(" +
            NUM_RE.pattern + r")\s*\(\s*rho\s*=\s*(" + NUM_RE.pattern + r")\s*\)", block_text, re.I
        )
        if lwl85m:
            vals = lwl85m.groups()
            row["LWL85m_dS"] = safe_float_token(vals[0])
            row["LWL85m_dN"] = safe_float_token(vals[1])
            row["LWL85m_omega"] = safe_float_token(vals[2])
            row["LWL85m_S"] = safe_float_token(vals[3])
            row["LWL85m_N"] = safe_float_token(vals[4])
            row["LWL85m_rho"] = safe_float_token(vals[5])

        lpb93 = re.search(
            r"LPB93\s*:\s*dS\s*=\s*(" + NUM_RE.pattern + r")\s*dN\s*=\s*(" + NUM_RE.pattern + r")\s*w\s*=\s*(" +
            NUM_RE.pattern + r")", block_text, re.I
        )
        if lpb93:
            vals = lpb93.groups()
            row["LPB93_dS"] = safe_float_token(vals[0])
            row["LPB93_dN"] = safe_float_token(vals[1])
            row["LPB93_omega"] = safe_float_token(vals[2])

    # attach global kappa if present (noninvasive)
    kappa_match = re.search(r'\bkappa\b\s*[:=]?\s*(' + NUM_RE.pattern + r')', content, re.I)
    kappa_global = safe_float_token(kappa_match.group(1)) if kappa_match else None
    if kappa_global is not None:
        for r in merged.values():
            r["kappa_global"] = kappa_global

    # Now augment each pair row with matrix values (if sibling files exist) and pairwise codon stats & QC flags
    # -------------------------------------------------------
    # DYNAMIC MATRIX HARVESTER (The "Maximalist" Upgrade)
    # -------------------------------------------------------
    base_dir = path.parent
    
    # We will store matrices as: matrix_store["2YN"]["dS"] = [vals...]
    matrix_store: Dict[str, Dict[str, List[float]]] = {}

    # Glob for all matrix-like files in the folder (2YN.dS, 2NG.dS, 2LWL.dS, etc.)
    # PAML naming convention is usually 2<METHOD>.<TYPE>
    matrix_files = list(base_dir.glob("2*.*"))
    
    for mfile in matrix_files:
        # Expecting filenames like '2YN.dS', '2NG.dN', '2LWL.t'
        parts = mfile.name.split('.')
        if len(parts) != 2: 
            continue
            
        method_prefix = parts[0]  # e.g., "2YN", "2NG", "2LWL"
        metric_type = parts[1]    # e.g., "dS", "dN", "t", "w" (some versions output w)
        
        # Only parse if it looks like a PAML matrix file (simple heuristic)
        if metric_type not in ("dS", "dN", "t", "w", "S", "N"):
            continue

        vals = parse_matrix_file_into_lower_tri_values(mfile)
        if vals:
            if method_prefix not in matrix_store:
                matrix_store[method_prefix] = {}
            matrix_store[method_prefix][metric_type] = vals

    # Helper to calculate matrix index (reused from your existing code)
    seqs = build_seq_token_list(content)
    nseq = len(seqs)
    pair_index_order: List[Tuple[str,str]] = []
    if nseq >= 2:
        for i in range(1, nseq):
            for j in range(0, i):
                pair_index_order.append((seqs[i], seqs[j]))
    pair_to_matrix_index = {}
    for idx, (a,b) in enumerate(pair_index_order):
        k = tuple(sorted([a,b]))
        pair_to_matrix_index[k] = idx
    final_rows = []
    for key, row in merged.items():
       # Inside the loop: for key, row in merged.items():
        idx = pair_to_matrix_index.get(key)
        
        # --- 1. Standard YN00 Matrices (Keep your existing columns) ---
        # We specifically look for '2YN' in our new dynamic store
        yn_vals = matrix_store.get("2YN", {})
        
        if idx is not None:
            row["matrix_dS"] = float(yn_vals["dS"][idx]) if (yn_vals.get("dS") and idx < len(yn_vals["dS"])) else ""
            row["matrix_dN"] = float(yn_vals["dN"][idx]) if (yn_vals.get("dN") and idx < len(yn_vals["dN"])) else ""
            row["matrix_t"]  = float(yn_vals["t"][idx])  if (yn_vals.get("t")  and idx < len(yn_vals["t"]))  else ""
        else:
            row["matrix_dS"] = ""
            row["matrix_dN"] = ""
            row["matrix_t"]  = ""

        # --- 2. NEW: Alternative Method Matrices (JSON Dump) ---
        # Captures 2NG (Nei-Gojobori), 2LWL, etc. for this specific pair
        alt_matrices_data = {}
        
        if idx is not None:
            for method, metrics in matrix_store.items():
                if method == "2YN": continue # Already handled
                
                alt_matrices_data[method] = {}
                for mtype, mvals in metrics.items():
                    if idx < len(mvals):
                        alt_matrices_data[method][mtype] = mvals[idx]
        
        row["alt_matrix_vals"] = json.dumps(alt_matrices_data) if alt_matrices_data else ""
        row["seq_labels"] = row.get("seq_labels", "")
        row["codon_counts"] = row.get("codon_counts", "")
        row["site_pattern_counts"] = row.get("site_pattern_counts", "")
        row["GC_pos1"] = row.get("GC_pos1", "")
        row["GC_pos2"] = row.get("GC_pos2", "")
        row["GC_pos3"] = row.get("GC_pos3", "")
        row["GC_total"] = row.get("GC_total", "")

        s1 = row.get("seq1")
        s2 = row.get("seq2")
        if before_seqmap and s1 and s2:
            idc, gaps, ambig, pct = compute_pair_codons_stats(before_seqmap, s1, s2)
            row["pair_n_identical_codons"] = idc
            row["pair_n_gapped_codons"] = gaps
            row["pair_n_ambiguous_codons"] = ambig
            row["pair_pct_identity"] = pct if pct is not None else ""
        else:
            row["pair_n_identical_codons"] = ""
            row["pair_n_gapped_codons"] = ""
            row["pair_n_ambiguous_codons"] = ""
            row["pair_pct_identity"] = ""

        ds_val = row.get("YN00_dS")
        dn_val = row.get("YN00_dN")
        dnse = row.get("YN00_dN_SE")
        dsse = row.get("YN00_dS_SE")
        omega_val = row.get("YN00_omega")

        row["dS_is_zero"] = (True if (ds_val == 0) else False) if ds_val is not None else False
        row["any_nan_inf"] = any(x is None for x in (ds_val, dn_val, dnse, dsse))
        try:
            row["omega_placeholder_flag"] = (True if (omega_val is not None and float(omega_val) >= 90.0) else False)
        except:
            row["omega_placeholder_flag"] = False

        # Ensure the new codon fields are present and stringified (JSON or empty)
        row["codon_usage_counts"] = row.get("codon_usage_counts", "") or codon_usage_counts_json or ""
        row["codon_pos_base_seq"] = row.get("codon_pos_base_seq", "") or codon_pos_base_seq_json or ""
        row["codon_pos_base_avg"] = row.get("codon_pos_base_avg", "") or codon_pos_base_avg_json or ""

        final_rows.append(row)

    return final_rows, f"✔ {gene_name}: Parsed {len(final_rows)} branches successfully"

# -------------------------
# Utility to format vector columns (unchanged)
# -------------------------
def _format_vector_from_row(row: Dict[str,Any], base_keys: List[str]) -> str:
    vals = []
    any_present = False
    for k in base_keys:
        v = row.get(k, None)
        if v is None:
            vals.append("NA")
        else:
            any_present = True
            vals.append(str(v))
    if not any_present:
        return ""
    return "[" + ", ".join(vals) + "]"

# ---------------- Public API ----------------
def run(input_path: str, output_path: str, *, verbose: bool = False) -> int:
    """
    Run yn00 harvester.

    input_path: a yn00.out file, a directory containing yn00.out, or a glob.
                If a directory, recursively finds yn00.out files within.
    output_path: CSV file path or directory. If directory, writes to
                 'yn00_filtered_stats.csv' inside it.
    Returns the number of rows written.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # Determine list of files
    if in_path.is_dir():
        # Recursively look for yn00.out files
        files = sorted(in_path.rglob("yn00.out"))
        if not files:
            # Fallback: try any .out file
            files = sorted(in_path.rglob("*.out"))
    elif in_path.is_file():
        files = [in_path]
    else:
        # Treat as glob
        files = sorted(Path(p) for p in glob.glob(str(in_path)))

    if not files:
        if verbose:
            print(f"[yn00] no input files found at {input_path}")
        return 0

    # Determine output file
    if out_path.is_dir():
        out_file = out_path / "yn00_filtered_stats.csv"
    else:
        out_file = out_path

    out_file.parent.mkdir(parents=True, exist_ok=True)

    rows_written = 0
    with out_file.open("w", newline="", encoding="utf-8") as outfh:
        writer = csv.DictWriter(outfh, fieldnames=FIELDNAMES, extrasaction='ignore')
        writer.writeheader()

        for fp in files:
            rows, msg = parse_yn00_file_rows(str(fp))
            if verbose and msg:
                print(f"[yn00] {msg}")

            for r in rows:
                # ... reuse the same row-formatting logic from write_streamed() ...
                outrow: Dict[str, Any] = {}
                for col in FIELDNAMES:
                    if col in ("L_i", "Ns_i", "Nv_i", "A_i", "B_i"):
                        if col in r and isinstance(r[col], str) and r[col].startswith('['):
                            outrow[col] = r[col]
                        else:
                            base_keys = {
                                "L_i": ["L_1", "L_2", "L_3"],
                                "Ns_i": ["Ns_1", "Ns_2", "Ns_3"],
                                "Nv_i": ["Nv_1", "Nv_2", "Nv_3"],
                                "A_i": ["A_1", "A_2", "A_3"],
                                "B_i": ["B_1", "B_2", "B_3"]
                            }[col]
                            outrow[col] = _format_vector_from_row(r, base_keys)
                        continue

                    if col in r:
                        outrow[col] = r.get(col, "")
                        continue
                    if col in BACKWARD_SCALAR_MAP:
                        orig = BACKWARD_SCALAR_MAP[col]
                        outrow[col] = r.get(orig, "")
                    else:
                        outrow[col] = r.get(col, "")

                outrow["codon_usage_counts"] = outrow.get("codon_usage_counts", "") or ""
                outrow["codon_pos_base_seq"] = outrow.get("codon_pos_base_seq", "") or ""
                outrow["codon_pos_base_avg"] = outrow.get("codon_pos_base_avg", "") or ""

                writer.writerow(outrow)
                rows_written += 1

    if verbose:
        print(f"[yn00] wrote {rows_written} rows to {out_file}")

    return rows_written


# ---------------- Standalone CLI (preserved for direct script use) ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="yn00_harvester", description="Harvest yn00 outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input yn00.out file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✔ yn00 stats written: {count} rows")
