#!/usr/bin/env python3
"""
13_basemlHarvester_with_pairwise.py (no alignment_frac columns)

- Run from v4_orthexplorer/ (CWD).
- Harvests baseml outputs under results/gene_centric/paml/*/baseml
- Produces ONE unified CSV:
    - baseml_unified_harvest_full.csv
      (branch-level rows with gene-level global params + pairwise kappa summaries)
"""
import re
import pandas as pd
import os
import glob
from Bio import Phylo, AlignIO
from io import StringIO
from pathlib import Path
import statistics
import json
import math

# ---------- helpers ----------
def float_or_none(s):
    try:
        return float(s)
    except Exception:
        return None

def strip_phylip_header_from_text(text):
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip() != ""]
    if lines and re.match(r"^\s*\d+\s+\d+", lines[0]):
        return "\n".join(lines[1:])
    return "\n".join(lines)

def find_tree_in_dir(bdir):
    p = Path(bdir)
    for name in ("tree_with_lengths_BS_label.nwk", "tree_with_lengths_BS.nwk", "tree.nwk"):
        f = p / name
        if f.exists():
            raw = f.read_text()
            return strip_phylip_header_from_text(raw)
    return None

# ---------- parse GTR and base frequencies ----------
def parse_gtr_and_pis(content):
    tail = content[-2000:] if len(content) > 2000 else content
    piA = piC = piG = piT = None
    mpi = re.search(r"Base frequencies\s*[:=]?\s*([\d\.\sEe+\-]+)", tail, re.I)
    if mpi:
        nums = re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?', mpi.group(1))
        if len(nums) >= 4:
            piA, piC, piG, piT = map(float, nums[:4])
    # GTR/REV five rates as in outputs
    rAC = rAG = rAT = rCG = rCT = None
    mgtr = re.search(r"Rate parameters\s*[:=]?\s*([\d\.\sEe+\-\s]+)", tail, re.I)
    if mgtr:
        nums = re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?', mgtr.group(1))
        if len(nums) >= 5:
            rAC, rAG, rAT, rCG, rCT = map(float, nums[:5])
    return (rAC, rAG, rAT, rCG, rCT), (piA, piC, piG, piT)

# ---------- parse global kappa and ts/tv ----------
def parse_global_kappa(content):
    tail = content[-2000:] if len(content) > 2000 else content
    patterns = [
        r'kappa\s*=\s*([-\d.eE]+)',
        r'kappa\s*\(ts/tv\)\s*=\s*([-\d.eE]+)',
        r'Rate matrix Q, Average Ts\/Tv\s*=\s*([-\d.eE]+)'
    ]
    for pat in patterns:
        m = re.search(pat, tail, re.I | re.S)
        if m:
            return float_or_none(m.group(1))
    return None

def parse_ts_tv_extracted(content):
    tail = content[-2000:] if len(content) > 2000 else content
    m = re.search(r'Average Ts\/Tv\s*=\s*([-\d.eE]+)', tail)
    if m:
        return float_or_none(m.group(1))
    m2 = re.search(r'(?i)(?:ts\/tv|transition\/transversion).*?([-\d.eE]+)', content)
    if m2:
        return float_or_none(m2.group(1))
    return None

# ---------- pairwise block extraction ----------
def extract_pairwise_block(content):
    lines = content.splitlines()
    start_idx = None
    for i, ln in enumerate(lines):
        if 'Distances' in ln and re.search(r'TN93', ln, re.I):
            start_idx = i + 1
            break
    if start_idx is None:
        for i, ln in enumerate(lines):
            if ln.strip().startswith('Distances:'):
                start_idx = i + 1
                break
    if start_idx is None:
        for i in range(len(lines)):
            if len(re.findall(r'\d+\.\d+', lines[i])) >= 2:
                start_idx = i
                break
    if start_idx is None:
        return []
    out_lines = []
    for ln in lines[start_idx:]:
        if re.match(r'^\s*(TREE\s+#|TREE\s|lnL\(|Detailed output|Parameters  in the rate matrix|TREE #|==>)', ln):
            break
        out_lines.append(ln.rstrip())
    while out_lines and out_lines[0].strip() == '':
        out_lines.pop(0)
    while out_lines and out_lines[-1].strip() == '':
        out_lines.pop()
    return out_lines

def parse_pairwise_matrix_lines(lines):
    parsed = []
    for ln in lines:
        if ln.strip() == '':
            continue
        m = re.match(r'^\s*(\S(?:.*?\S)?)\s{2,}(.+)$', ln)
        if m:
            label = m.group(1).strip()
            rest = m.group(2).strip()
        else:
            m2 = re.match(r'^\s*([^\d:()]+?)\s+([0-9\.\-].*)$', ln)
            if m2:
                label = m2.group(1).strip()
                rest = m2.group(2).strip()
            else:
                if parsed:
                    prev_label, prev_rest = parsed[-1]
                    rest = (prev_rest + ' ' + ln.strip()).strip()
                    parsed[-1] = (prev_label, rest)
                continue
        parsed.append((label, rest))

    labels = [p[0] for p in parsed]
    full_pairs = {}
    for i, (lab, rest) in enumerate(parsed):
        pair_pattern = re.findall(r'([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)(?:\s*\(\s*([-\d\.]+)\s*\))?', rest)
        if pair_pattern:
            for j_idx, (dist_s, kappa_s) in enumerate(pair_pattern):
                j = j_idx
                if j >= len(labels):
                    break
                seq_i = labels[i]
                seq_j = labels[j]
                d = float_or_none(dist_s)
                k = float_or_none(kappa_s) if kappa_s else None
                full_pairs[(seq_i, seq_j)] = {'distance': d, 'kappa': k}
                if (seq_j, seq_i) not in full_pairs:
                    full_pairs[(seq_j, seq_i)] = {'distance': d, 'kappa': k}
        else:
            nums = re.findall(r'([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)', rest)
            for j_idx, n in enumerate(nums):
                j = j_idx
                if j >= len(labels):
                    break
                seq_i = labels[i]
                seq_j = labels[j]
                d = float_or_none(n)
                full_pairs[(seq_i, seq_j)] = {'distance': d, 'kappa': None}
                if (seq_j, seq_i) not in full_pairs:
                    full_pairs[(seq_j, seq_i)] = {'distance': d, 'kappa': None}
    return labels, full_pairs

# ---------- additional parsers ----------
def extract_constant_sites(content):
    m = re.search(r'#\s*constant sites\s*[:=]?\s*(\d+)\s*\(\s*([\d\.]+)%\s*\)', content, re.I)
    if not m:
        m = re.search(r'constant sites[:\s]+(\d+)\s*\(?([\d\.]+)%\)?', content, re.I)
    if m:
        return int(m.group(1)), float(m.group(2))
    return None, None

def extract_homogeneity(content):
    m = re.search(r'Homogeneity statistic[:\s]*X2\s*=\s*([-\d.eE]+)\s*G\s*=\s*([-\d.eE]+)', content, re.I)
    if m:
        return float_or_none(m.group(1)), float_or_none(m.group(2))
    return None, None

def extract_mp_score(content):
    m = re.search(r'MP score[:\s]*([-\d\.]+)', content)
    return float_or_none(m.group(1)) if m else None

def extract_n_sites_from_content(content):
    lines = content.splitlines()
    if lines:
        m = re.match(r'^\s*(\d+)\s+(\d+)\s*$', lines[0])
        if m:
            return int(m.group(2)), int(m.group(1))
    m2 = re.search(r'(?:Number of sites|sites)[:\s]+([0-9]+)', content, re.I)
    if m2:
        return int(m2.group(1)), None
    return None, None

def extract_rate_matrix_Q(content):
    m = re.search(r'Rate matrix Q[^\n]*\n((?:\s*[-\d\.\sE+]+\n){4})', content, re.I)
    if not m:
        return None
    block = m.group(1).strip().splitlines()
    rows = []
    for r in block:
        nums = re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?', r)
        rows.append([float_or_none(x) for x in nums])
    if len(rows) == 4:
        return json.dumps(rows)
    return None

def parse_lnL(content):
    m = re.search(r'lnL\([^\)]*\)\s*[:=]?\s*([-\d\.eE]+)', content)
    if m:
        return float_or_none(m.group(1))
    m2 = re.search(r'Negative\s*log\s*likelihood\s*=\s*([-\d\.eE]+)', content, re.I)
    if m2:
        return float_or_none(m2.group(1))
    m3 = re.findall(r'-\d{3,}\.\d+', content)
    if m3:
        return float_or_none(m3[-1])
    return None

def compute_info_criteria(lnL, k, n):
    if lnL is None or k is None:
        return None, None, None
    AIC = 2.0 * k - 2.0 * lnL
    AICc = None
    if n is not None and (n - k - 1) > 0:
        AICc = AIC + (2.0 * k * (k + 1.0)) / (n - k - 1.0)
    BIC = k * math.log(n) - 2.0 * lnL if n is not None and n > 0 else None
    return AIC, AICc, BIC

# ---------- compute tree metrics ----------
def compute_tree_metrics_from_newick(newick_text):
    if not newick_text:
        return {}, None
    try:
        tree = Phylo.read(StringIO(strip_phylip_header_from_text(newick_text)), "newick")
    except Exception as e:
        print(f"[WARN] compute_tree_metrics: phylo read failed: {e}")
        return {}, None
    counter = 0
    for cl in tree.find_clades(order='preorder'):
        if cl.is_terminal():
            continue
        if cl.name is None:
            cl.name = f"Node{counter}"
            counter += 1
    parent_map = {}
    for parent in tree.find_clades(order='preorder'):
        for child in parent.clades:
            parent_map[child] = parent
    root = tree.root
    root_to_dist = {}
    def fill_dist(node, acc=0.0):
        root_to_dist[node] = acc
        for c in node.clades:
            bl = c.branch_length if c.branch_length is not None else 0.0
            fill_dist(c, acc + bl)
    fill_dist(root, 0.0)
    metrics = {}
    for cl in tree.find_clades(order='preorder'):
        is_internal = not cl.is_terminal()
        n_desc = len(cl.get_terminals()) if is_internal else 1
        parent = parent_map.get(cl)
        parent_name = parent.name if parent is not None else None
        rt = root_to_dist.get(cl, None)
        support = None
        try:
            if hasattr(cl, 'confidence'):
                support = float_or_none(cl.confidence)
        except Exception:
            support = None
        metrics[cl.name] = {
            'parent': parent_name,
            'is_internal': is_internal,
            'n_descendants': n_desc,
            'root_to_tip': rt,
            'support': support
        }
    n_taxa = len(tree.get_terminals())
    total_tree_length = sum((c.branch_length or 0.0) for c in tree.find_clades())
    return metrics, (n_taxa, total_tree_length)

# ---------- harvesting ----------
def harvest_baseml_branch_and_pairs(gene_id, baseml_path):
    if not os.path.exists(baseml_path):
        return [], [], None

    with open(baseml_path, 'r') as fh:
        content = fh.read()

    rates, pis = parse_gtr_and_pis(content)
    rAC, rAG, rAT, rCG, rCT = rates
    piA, piC, piG, piT = pis
    global_kappa = parse_global_kappa(content)
    ts_tv_global = parse_ts_tv_extracted(content)
    constant_sites, constant_pct = extract_constant_sites(content)
    hom_X2, hom_G = extract_homogeneity(content)
    mp_score = extract_mp_score(content)
    n_sites_guess, n_taxa_guess = extract_n_sites_from_content(content)
    ntime = np_param = None
    m_time = re.search(r'lnL\([^\)]*ntime[:\s]*([0-9]+)[^\)]*np[:\s]*([0-9]+)\)', content)
    if m_time:
        try:
            ntime = int(m_time.group(1))
            np_param = int(m_time.group(2))
        except Exception:
            ntime = None
            np_param = None
    Q_matrix = extract_rate_matrix_Q(content)
    baseml_dir = os.path.dirname(baseml_path)
    newick = find_tree_in_dir(baseml_dir)
    tree_metrics_map, tree_meta = compute_tree_metrics_from_newick(newick)
    n_taxa_from_tree = tree_meta[0] if tree_meta else None
    total_tree_length_from_tree = tree_meta[1] if tree_meta else None
    branch_rows = []
    node_list = []
    if newick:
        try:
            t = Phylo.read(StringIO(strip_phylip_header_from_text(newick)), "newick")
            node_list = [c for c in t.find_clades(order='preorder')]
            bls = [c.branch_length for c in node_list if c.branch_length is not None]
            tree_length = sum(bls) if bls else (total_tree_length_from_tree if total_tree_length_from_tree is not None else None)
        except Exception as e:
            print(f"[WARN] tree parse {gene_id}: {e}")
            node_list = []
            tree_length = total_tree_length_from_tree
    else:
        tree_length = total_tree_length_from_tree
    lnL = parse_lnL(content)
    if node_list:
        for idx, clade in enumerate(node_list):
            if clade.branch_length is None:
                continue
            name = clade.name if clade.name else f"Node{idx}"
            tmetrics = tree_metrics_map.get(name, {})
            row = {
                'Gene_ID': gene_id,
                'Branch_ID': name,
                'Branch_Length_t': clade.branch_length,
                'lnL': lnL,
                'Alpha_Fixed': 0.5,
                'GTR_a': rAC,
                'GTR_b': rAG,
                'GTR_c': rAT,
                'GTR_d': rCG,
                'GTR_e': rCT,
                'piA': piA,
                'piC': piC,
                'piG': piG,
                'piT': piT,
                'kappa': global_kappa,
                'kappa_source': "global_model" if global_kappa is not None else "missing",
                'ts_tv': ts_tv_global,
                'ts_tv_source': "global_model" if ts_tv_global is not None else "missing",
                'tree_length': tree_length,
                'parent': tmetrics.get('parent'),
                'is_internal': tmetrics.get('is_internal'),
                'n_descendants': tmetrics.get('n_descendants'),
                'root_to_tip': tmetrics.get('root_to_tip'),
                'branch_support': tmetrics.get('support')
            }
            branch_rows.append(row)
    pairwise_rows = []
    pair_lines = extract_pairwise_block(content)
    labels = []
    full_pairs = {}
    if pair_lines:
        labels, full_pairs = parse_pairwise_matrix_lines(pair_lines)
        for (s1, s2), v in full_pairs.items():
            pairwise_rows.append({
                'Gene_ID': gene_id,
                'Seq1': s1,
                'Seq2': s2,
                'distance': v.get('distance'),
                'kappa': v.get('kappa')
            })
    AIC, AICc, BIC = compute_info_criteria(lnL, np_param, n_sites_guess)
    gene_summary = {
        'Gene_ID': gene_id,
        'n_sequences': len(labels) if labels else (n_taxa_from_tree if n_taxa_from_tree is not None else None),
        'n_pairs': len(pairwise_rows),
        'pairwise_kappa_mean': None,
        'pairwise_kappa_median': None,
        'pairwise_kappa_sd': None,
        'pairwise_kappa_min': None,
        'pairwise_kappa_max': None,
        'pairwise_distance_mean': None,
        'pairwise_distance_min': None,
        'pairwise_distance_max': None,
        'global_kappa': global_kappa,
        'ts_tv': ts_tv_global,
        'tree_length': tree_length,
        'n_taxa': n_taxa_from_tree if n_taxa_from_tree is not None else n_taxa_guess,
        'n_sites': n_sites_guess,
        'constant_sites': constant_sites,
        'constant_pct': constant_pct,
        'hom_X2': hom_X2,
        'hom_G': hom_G,
        'ntime': ntime,
        'np_param': np_param,
        'mp_score': mp_score,
        'Q_matrix': Q_matrix,
        'AIC': AIC,
        'AICc': AICc,
        'BIC': BIC
    }
    kappas = [r['kappa'] for r in pairwise_rows if r['kappa'] is not None]
    dists = [r['distance'] for r in pairwise_rows if r['distance'] is not None]
    if kappas:
        gene_summary['pairwise_kappa_mean'] = statistics.mean(kappas)
        gene_summary['pairwise_kappa_median'] = statistics.median(kappas)
        gene_summary['pairwise_kappa_sd'] = statistics.pstdev(kappas) if len(kappas) > 1 else 0.0
        gene_summary['pairwise_kappa_min'] = min(kappas)
        gene_summary['pairwise_kappa_max'] = max(kappas)
    if dists:
        gene_summary['pairwise_distance_mean'] = statistics.mean(dists)
        gene_summary['pairwise_distance_min'] = min(dists)
        gene_summary['pairwise_distance_max'] = max(dists)
    return branch_rows, pairwise_rows, gene_summary

# ---------------- Public API ----------------
def run(input_path: str, output_path: str, *, verbose: bool = False) -> int:
    """
    Run baseml harvester.

    input_path: a baseml.out file, a directory containing baseml.out files,
                or a glob pattern. If a directory, recursively finds baseml.out.
    output_path: CSV file path or directory. If directory, writes to
                 'baseml_filtered_stats.csv' inside it.
    Returns the number of rows written.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # Determine list of baseml.out files
    if in_path.is_dir():
        files = sorted(in_path.rglob("baseml.out"))
    elif in_path.is_file():
        files = [in_path]
    else:
        files = sorted(Path(p) for p in glob.glob(str(in_path)))

    if not files:
        if verbose:
            print(f"[baseml] no input files found at {input_path}")
        return 0

    # Determine output file
    if out_path.is_dir():
        out_file = out_path / "baseml_filtered_stats.csv"
    else:
        out_file = out_path

    out_file.parent.mkdir(parents=True, exist_ok=True)

    all_branch_rows = []
    all_pairwise_rows = []
    all_gene_summaries = []

    for f_path in files:
        # Gene name: parent directory of the baseml/ subdirectory
        # Conventional layout: <gene>/baseml/baseml.out
        try:
            gene_name = Path(f_path).parents[1].name
        except IndexError:
            gene_name = Path(f_path).stem

        try:
            brows, prows, gsum = harvest_baseml_branch_and_pairs(gene_name, str(f_path))
            if brows:
                all_branch_rows.extend(brows)
            if prows:
                all_pairwise_rows.extend(prows)
            if gsum:
                all_gene_summaries.append(gsum)
            if verbose:
                print(f"[baseml] processed {gene_name} (branches: {len(brows)}, pairs: {len(prows)})")
        except Exception as e:
            if verbose:
                print(f"[baseml] ERROR {gene_name}: {e}")

    if not all_branch_rows:
        if verbose:
            print("[baseml] no branch-level rows produced")
        return 0

    # Build the merged DataFrame (same logic as the original __main__ block)
    df_br = pd.DataFrame(all_branch_rows)

    if all_pairwise_rows:
        df_pairs = pd.DataFrame(all_pairwise_rows)
        df_pairs_unique = df_pairs[df_pairs.apply(lambda r: (r['Seq1'] <= r['Seq2']), axis=1)].copy()
        pw_summary = df_pairs_unique.groupby("Gene_ID")["kappa"].agg(
            pairwise_kappa_mean="mean", pairwise_kappa_min="min", pairwise_kappa_max="max"
        ).reset_index()
        pw_dist_summary = df_pairs_unique.groupby("Gene_ID")["distance"].agg(
            pairwise_distance_mean="mean", pairwise_distance_min="min", pairwise_distance_max="max"
        ).reset_index()
        pw_summary = pw_summary.merge(pw_dist_summary, on="Gene_ID", how="outer")
    else:
        pw_summary = pd.DataFrame(columns=[
            "Gene_ID", "pairwise_kappa_mean", "pairwise_kappa_min", "pairwise_kappa_max",
            "pairwise_distance_mean", "pairwise_distance_min", "pairwise_distance_max"
        ])

    df_gsum = pd.DataFrame(all_gene_summaries) if all_gene_summaries else pd.DataFrame(columns=[
        'Gene_ID', 'global_kappa', 'ts_tv', 'tree_length', 'n_taxa', 'n_sites',
        'constant_sites', 'constant_pct', 'hom_X2', 'hom_G',
        'ntime', 'np_param', 'mp_score', 'Q_matrix', 'AIC', 'AICc', 'BIC',
        'pairwise_distance_mean', 'pairwise_distance_min', 'pairwise_distance_max',
        'pairwise_kappa_mean', 'pairwise_kappa_min', 'pairwise_kappa_max'
    ])

    df_merged = df_br.merge(pw_summary, on="Gene_ID", how="left")

    if not df_gsum.empty:
        gsum_cols = [
            'Gene_ID', 'global_kappa', 'ts_tv', 'tree_length', 'n_taxa', 'n_sites',
            'constant_sites', 'constant_pct', 'hom_X2', 'hom_G',
            'ntime', 'np_param', 'mp_score', 'Q_matrix', 'AIC', 'AICc', 'BIC',
            'pairwise_distance_mean', 'pairwise_distance_min', 'pairwise_distance_max'
        ]
        gsum_small = df_gsum[[c for c in gsum_cols if c in df_gsum.columns]].copy()
        df_merged = df_merged.merge(gsum_small, on='Gene_ID', how='left', suffixes=('', '_gs'))

        df_merged['kappa'] = df_merged.apply(
            lambda r: r['kappa'] if pd.notnull(r['kappa']) else (
                r['global_kappa'] if pd.notnull(r.get('global_kappa')) else None
            ),
            axis=1
        )
        df_merged['kappa_source'] = df_merged.apply(
            lambda r: r.get('kappa_source') if (
                r.get('kappa_source') and r.get('kappa_source') != 'missing'
                and pd.notnull(r.get('kappa'))
            ) else (
                'global_model' if pd.notnull(r.get('global_kappa'))
                else (r.get('kappa_source') if 'kappa_source' in r else 'missing')
            ),
            axis=1
        )
        df_merged['ts_tv'] = df_merged.apply(
            lambda r: r['ts_tv'] if pd.notnull(r['ts_tv']) else (
                r['ts_tv_gs'] if pd.notnull(r.get('ts_tv_gs')) else None
            ),
            axis=1
        )
        df_merged['ts_tv_source'] = df_merged.apply(
            lambda r: r.get('ts_tv_source') if (
                r.get('ts_tv_source') and r.get('ts_tv_source') != 'missing'
                and pd.notnull(r.get('ts_tv'))
            ) else (
                'global_model' if pd.notnull(r.get('ts_tv_gs'))
                else (r.get('ts_tv_source') if 'ts_tv_source' in r else 'missing')
            ),
            axis=1
        )
        df_merged = df_merged.drop(
            columns=[c for c in ['global_kappa', 'ts_tv_gs', 'tree_length_gs'] if c in df_merged.columns],
            errors='ignore'
        )

    final_cols = [
        'Gene_ID', 'Branch_ID', 'Branch_Length_t', 'lnL', 'Alpha_Fixed',
        'GTR_a', 'GTR_b', 'GTR_c', 'GTR_d', 'GTR_e',
        'piA', 'piC', 'piG', 'piT',
        'kappa', 'kappa_source', 'ts_tv', 'ts_tv_source',
        'pairwise_kappa_mean', 'pairwise_kappa_min', 'pairwise_kappa_max',
        'pairwise_distance_mean', 'pairwise_distance_min', 'pairwise_distance_max',
        'tree_length',
        'parent', 'is_internal', 'n_descendants', 'root_to_tip', 'branch_support',
        'n_taxa', 'n_sites', 'constant_sites', 'constant_pct', 'hom_X2', 'hom_G',
        'ntime', 'np_param', 'mp_score', 'Q_matrix',
        'AIC', 'AICc', 'BIC'
    ]
    final_existing = [c for c in final_cols if c in df_merged.columns]
    df_final = df_merged[final_existing].copy()

    df_final.to_csv(out_file, index=False)

    if verbose:
        print(f"[baseml] wrote {len(df_final)} rows to {out_file}")

    return len(df_final)


# ---------------- Standalone CLI (preserved for direct script use) ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="baseml_harvester", description="Harvest baseml outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input baseml.out file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✔ baseml stats written: {count} rows")
