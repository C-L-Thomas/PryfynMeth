#!/usr/bin/env python3
"""
CpG Context Counter (+ per-site summary, optional detailed file,
and Fisher tests incl. left/right aggregates with enrichment)

Usage
-----
# Genome-wide counts only
python cpg_context_counter.py genome.fa -o out.tsv

# Counts + per-site SUMMARY appended below + Fisher tests
python cpg_context_counter.py genome.fa --sites sites.tsv -o out.tsv

# Also write per-site DETAILED rows to a separate file
python cpg_context_counter.py genome.fa --sites sites.tsv --sites-detailed-out sites_detailed.tsv -o out.tsv

# Choose Fisher results output name
python cpg_context_counter.py genome.fa --sites sites.tsv --fishers-out my_fisher.tsv -o out.tsv
"""

from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from itertools import product
from typing import Dict, Generator, Iterable, Tuple, List, Set
import sys
import csv
import math

DNA = set("ACGT")

# ---------------------------- FASTA helpers ---------------------------- #

def parse_fasta(path: str) -> Generator[Tuple[str, str], None, None]:
    """Yield (header, sequence) for each record in a FASTA file. Uppercases sequence."""
    header = None
    seq_chunks = []
    with open(path, 'r') as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_chunks).upper()
                header = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(seq_chunks).upper()

def revcomp(seq: str) -> str:
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(comp)[::-1]

def all_context_keys(flank: int, collapse_rc: bool) -> Iterable[str]:
    alphabet = ['A', 'C', 'G', 'T']
    if collapse_rc:
        seen = set()
        for left in product(alphabet, repeat=flank):
            for right in product(alphabet, repeat=flank):
                ctx = ''.join(left) + 'CG' + ''.join(right)
                rc = revcomp(ctx)
                key = min(ctx, rc)
                if key not in seen:
                    seen.add(key)
                    yield key
    else:
        for left in product(alphabet, repeat=flank):
            for right in product(alphabet, repeat=flank):
                yield ''.join(left) + 'CG' + ''.join(right)

# ---------------------------- Genome-wide counts ---------------------------- #

def scan_sequence(seq: str, flank: int, include_ambig: bool, collapse_rc: bool):
    """
    Return:
      context_counts: Counter
      stats: dict
      side_counts: {'left': Counter({'A':..,'C':..,'G':..,'T':..}),
                    'right': Counter({...})}
    Side tallies use the immediate neighbor bases (one base) in plus orientation.
    """
    context_counts: Counter = Counter()
    stats = defaultdict(int)
    side_counts = {'left': Counter(), 'right': Counter()}
    n = len(seq)
    i = 0
    while True:
        i = seq.find('CG', i)
        if i == -1:
            break
        stats['total_cpg_sites'] += 1
        left_start = i - flank
        right_end = i + 2 + flank
        if left_start < 0 or right_end > n:
            stats['skipped_edge'] += 1
            i += 1
            continue
        left = seq[left_start:i]
        right = seq[i+2:right_end]
        if not include_ambig and (set(left) - DNA or set(right) - DNA):
            stats['skipped_ambiguous'] += 1
            i += 1
            continue
        ctx = left + 'CG' + right
        if collapse_rc:
            ctx = min(ctx, revcomp(ctx))
        context_counts[ctx] += 1
        stats['kept'] += 1
        if flank >= 1:
            lb = left[0] if left else None
            rb = right[0] if right else None
            if lb in DNA:
                side_counts['left'][lb] += 1
            if rb in DNA:
                side_counts['right'][rb] += 1
        i += 1
    return context_counts, stats, side_counts

# ---------------------------- Sites helpers ---------------------------- #

def sniff_delimiter(first_line: str) -> str | None:
    if '\t' in first_line:
        return '\t'
    if ',' in first_line:
        return ','
    return None

def load_sites(path: str) -> Tuple[List[dict], Set[str]]:
    """Load a sites table. Requires columns: chr, pos, strand. Returns (rows, set_of_chromosomes)."""
    rows: List[dict] = []
    chroms: Set[str] = set()
    with open(path, 'r', newline='') as f:
        peek = f.readline()
        if not peek:
            return rows, chroms
        delim = sniff_delimiter(peek)
        f.seek(0)
        if delim:
            reader = csv.DictReader(f, delimiter=delim)
            for r in reader:
                if not r:
                    continue
                if 'chr' not in r or 'pos' not in r or 'strand' not in r:
                    raise ValueError("Sites file must contain headers: chr, pos, strand")
                r['pos'] = int(r['pos'])
                rows.append(r)
                chroms.add(r['chr'])
        else:
            header = peek.strip().split()
            if {'chr','pos','strand'} - set(header):
                raise ValueError("Sites file must contain headers: chr, pos, strand")
            idx = {h:i for i,h in enumerate(header)}
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                r = {h: parts[idx[h]] for h in header if idx[h] < len(parts)}
                r['pos'] = int(r['pos'])
                rows.append(r)
                chroms.add(r['chr'])
    return rows, chroms

def extract_context(seq: str, pos_1based: int, strand: str, flank: int) -> Tuple[str, str, str, str]:
    """
    '+' : CpG at [pos, pos+1] (1-based)
    '-' : CpG at [pos-1, pos] on '+' reference, then reverse-complement window so CG is centered on '-' orientation.
    Returns (context, left, right, status). left/right are immediate flanks on the returned orientation.
    """
    n = len(seq)
    if strand == '+':
        s = pos_1based - 1
        left_start = s - flank
        right_end = s + 2 + flank
        if left_start < 0 or right_end > n:
            return 'NA', 'NA', 'NA', 'out_of_bounds'
        context = seq[left_start:right_end]
        left = seq[left_start:s]
        right = seq[s+2:right_end]
        status = 'ok' if seq[s:s+2] == 'CG' else 'ref_not_CG'
        return context, left, right, status
    else:
        c = pos_1based - 1        # index of minus-strand C on + reference (this is a G on +)
        s = c - 1                 # start of CG on + reference should be at c-1
        left_start = s - flank
        right_end = s + 2 + flank
        if left_start < 0 or right_end > n:
            return 'NA', 'NA', 'NA', 'out_of_bounds'
        window_plus = seq[left_start:right_end]
        status = 'ok' if seq[s:s+2] == 'CG' else 'ref_not_CG'
        context_minus = revcomp(window_plus)
        left = context_minus[:flank]
        right = context_minus[-flank:] if flank > 0 else ''
        return context_minus, left, right, status

# ---------------------------- Per-site summary & detailed ---------------------------- #

def summarize_sites_contexts(site_rows: List[dict],
                             chrom_seqs: Dict[str, str],
                             flank: int,
                             collapse_rc: bool) -> Tuple[Counter, Dict[str,int], Dict[str,Counter]]:
    """
    Returns:
      counts (per-context),
      stats,
      side_counts {'left':Counter, 'right':Counter} using immediate neighbors on oriented context.
    """
    counts = Counter()
    stats = defaultdict(int)
    side_counts = {'left': Counter(), 'right': Counter()}
    for r in site_rows:
        chrom = r['chr']
        pos = int(r['pos'])
        strand = r['strand']
        seq = chrom_seqs.get(chrom)
        stats['total_sites'] += 1
        if seq is None:
            stats['chrom_missing'] += 1
            continue
        ctx, left, right, status = extract_context(seq, pos, strand, flank)
        if status != 'ok':
            stats[status] += 1
            continue
        if collapse_rc:
            ctx = min(ctx, revcomp(ctx))
        counts[ctx] += 1
        stats['kept_ok'] += 1
        if flank >= 1:
            lb = left[0] if left else None
            rb = right[0] if right else None
            if lb in DNA:
                side_counts['left'][lb] += 1
            if rb in DNA:
                side_counts['right'][rb] += 1
    return counts, stats, side_counts

def write_sites_detailed(path: str,
                         site_rows: List[dict],
                         chrom_seqs: Dict[str, str],
                         flank: int) -> None:
    context_len = 2 + 2*flank
    with open(path, 'w') as out:
        orig_cols = list(site_rows[0].keys()) if site_rows else ['chr','pos','strand']
        extra = [f'context_{context_len}bp','left','right','status']
        out.write('\t'.join(orig_cols + extra) + '\n')
        for r in site_rows:
            chrom = r['chr']
            pos = int(r['pos'])
            strand = r['strand']
            seq = chrom_seqs.get(chrom)
            if seq is None:
                context, left, right, status = ('NA','NA','NA','chrom_missing')
            else:
                context, left, right, status = extract_context(seq, pos, strand, flank)
            values = [str(r.get(col, '')) for col in orig_cols] + [context, left, right, status]
            out.write('\t'.join(values) + '\n')

# ---------------------------- Fisher's exact test (pure Python) ---------------------------- #

def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float('-inf')
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def _hypergeom_logpmf(a: int, r1: int, c1: int, n: int) -> float:
    return _log_comb(c1, a) + _log_comb(n - c1, r1 - a) - _log_comb(n, r1)

def _fisher_pvalues(a: int, b: int, c: int, d: int) -> Tuple[float, float, float, float]:
    """Return (odds_ratio_raw, p_less, p_greater, p_two_sided)."""
    or_raw = (a * d) / (b * c) if (b > 0 and c > 0) else (
        float('inf') if (a > 0 and d > 0 and (b == 0 or c == 0)) else 0.0
    )
    r1 = a + b
    r2 = c + d
    c1 = a + c
    c2 = b + d
    n = r1 + r2
    a_min = max(0, r1 - c2)
    a_max = min(r1, c1)
    logp_obs = _hypergeom_logpmf(a, r1, c1, n)
    p_less = 0.0
    for aa in range(a_min, a + 1):
        p_less += math.exp(_hypergeom_logpmf(aa, r1, c1, n))
    p_greater = 0.0
    for aa in range(a, a_max + 1):
        p_greater += math.exp(_hypergeom_logpmf(aa, r1, c1, n))
    p_two = 0.0
    for aa in range(a_min, a_max + 1):
        lp = _hypergeom_logpmf(aa, r1, c1, n)
        if lp <= logp_obs + 1e-15:
            p_two += math.exp(lp)
    # Clamp
    p_less = min(max(p_less, 0.0), 1.0)
    p_greater = min(max(p_greater, 0.0), 1.0)
    p_two = min(max(p_two, 0.0), 1.0)
    return or_raw, p_less, p_greater, p_two

def _benjamini_hochberg(pvals: List[float]) -> List[float]:
    """Return BH-FDR q-values in the same order as pvals."""
    m = len(pvals)
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    q = [0.0] * m
    prev = 1.0
    for rank, (i, p) in enumerate(indexed, start=1):
        val = p * m / rank
        prev = min(prev, val)
        q[i] = prev
    return q

def _safe_enrichment(a: int, fg_total: int, c: int, bg_total: int, pseudo: float = 0.5) -> Tuple[float, float]:
    """Risk ratio with light smoothing (Laplace add-0.5) to avoid divide-by-zero."""
    if fg_total <= 0 or bg_total <= 0:
        return float('nan'), float('nan')
    rate_fg = (a + pseudo) / (fg_total + 1.0)
    rate_bg = (c + pseudo) / (bg_total + 1.0)
    enr = rate_fg / rate_bg
    try:
        log2_enr = math.log2(enr)
    except ValueError:
        log2_enr = float('nan')
    return enr, log2_enr

def _write_fisher_section(out, rows, pvals_two, title, include_left: bool, include_right: bool):
    out.write(f"\n## {title}\n")
    headers = ['feature']
    if include_left:
        headers.append('left_base')
    if include_right:
        headers.append('right_base')
    headers += [
        'methylated_count','methylated_total',
        'background_count','background_total',
        'enrichment','log2_enrichment',
        'odds_ratio','odds_ratio_haldane',
        'p_less','p_greater','p_two_sided','fdr_two_sided'
    ]
    out.write('\t'.join(headers) + '\n')
    qvals = _benjamini_hochberg(pvals_two) if pvals_two else []
    for (row, q) in zip(rows, qvals):
        fields = [str(row['label'])]
        if include_left:
            fields.append(str(row['left']))
        if include_right:
            fields.append(str(row['right']))
        fields += [
            str(row['fg_count']), str(row['fg_total']),
            str(row['bg_count']), str(row['bg_total']),
            f"{row['enrichment']:.6g}", f"{row['log2_enrichment']:.6g}",
            f"{row['or_raw']:.6g}", f"{row['or_hal']:.6g}",
            f"{row['p_less']:.6g}", f"{row['p_greater']:.6g}",
            f"{row['p_two']:.6g}", f"{q:.6g}"
        ]
        out.write('\t'.join(fields) + '\n')

def run_fishers(fg_counts: Counter,
                bg_counts: Counter,
                flank: int,
                collapse_rc: bool,
                fishers_out: str,
                fg_sides: Dict[str, Counter] | None = None,
                bg_sides: Dict[str, Counter] | None = None) -> None:
    keys = list(all_context_keys(flank, collapse_rc))
    for k in keys:
        fg_counts.setdefault(k, 0)
        bg_counts.setdefault(k, 0)
    fg_total = sum(fg_counts.values())
    bg_total = sum(bg_counts.values())

    with open(fishers_out, 'w') as out:
        # --- Per-context section (keeps left & right columns) ---
        rows = []
        p2 = []
        for k in keys:
            a = fg_counts[k]; b = fg_total - a
            c = bg_counts[k]; d = bg_total - c
            or_raw, p_less, p_greater, p_two = _fisher_pvalues(a, b, c, d)
            or_hal = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            enr, log2_enr = _safe_enrichment(a, fg_total, c, bg_total)
            rows.append({
                'label': k, 'left': k[:flank], 'right': (k[-flank:] if flank>0 else ''),
                'fg_count': a, 'fg_total': fg_total,
                'bg_count': c, 'bg_total': bg_total,
                'enrichment': enr, 'log2_enrichment': log2_enr,
                'or_raw': or_raw, 'or_hal': or_hal,
                'p_less': p_less, 'p_greater': p_greater, 'p_two': p_two
            })
            p2.append(p_two)
        _write_fisher_section(out, rows, p2, 'Per-context Fisher tests', include_left=True, include_right=True)

        # --- Left/right aggregates (drop the opposite base column) ---
        if flank >= 1 and fg_sides is not None and bg_sides is not None:
            # LEFT: no right_base column
            rows = []; p2 = []
            fgL = fg_sides.get('left', Counter()); bgL = bg_sides.get('left', Counter())
            fgL_total = sum(fgL.values());        bgL_total = sum(bgL.values())
            for base in ['A','C','G','T']:
                a = fgL.get(base, 0); b = fgL_total - a
                c = bgL.get(base, 0); d = bgL_total - c
                or_raw, p_less, p_greater, p_two = _fisher_pvalues(a, b, c, d)
                or_hal = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                enr, log2_enr = _safe_enrichment(a, fgL_total, c, bgL_total)
                rows.append({
                    'label': f'left_{base}', 'left': base, 'right': '',
                    'fg_count': a, 'fg_total': fgL_total,
                    'bg_count': c, 'bg_total': bgL_total,
                    'enrichment': enr, 'log2_enrichment': log2_enr,
                    'or_raw': or_raw, 'or_hal': or_hal,
                    'p_less': p_less, 'p_greater': p_greater, 'p_two': p_two
                })
                p2.append(p_two)
            _write_fisher_section(out, rows, p2, 'Left-base aggregate Fisher tests', include_left=True, include_right=False)

            # RIGHT: no left_base column
            rows = []; p2 = []
            fgR = fg_sides.get('right', Counter()); bgR = bg_sides.get('right', Counter())
            fgR_total = sum(fgR.values());         bgR_total = sum(bgR.values())
            for base in ['A','C','G','T']:
                a = fgR.get(base, 0); b = fgR_total - a
                c = bgR.get(base, 0); d = bgR_total - c
                or_raw, p_less, p_greater, p_two = _fisher_pvalues(a, b, c, d)
                or_hal = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                enr, log2_enr = _safe_enrichment(a, fgR_total, c, bgR_total)
                rows.append({
                    'label': f'right_{base}', 'left': '', 'right': base,
                    'fg_count': a, 'fg_total': fgR_total,
                    'bg_count': c, 'bg_total': bgR_total,
                    'enrichment': enr, 'log2_enrichment': log2_enr,
                    'or_raw': or_raw, 'or_hal': or_hal,
                    'p_less': p_less, 'p_greater': p_greater, 'p_two': p_two
                })
                p2.append(p_two)
            _write_fisher_section(out, rows, p2, 'Right-base aggregate Fisher tests', include_left=False, include_right=True)
        else:
            out.write("\n## Left/right aggregate Fisher tests\n")
            out.write("# Skipped (requires flank>=1 and valid side tallies)\n")

# ---------------------------- Output (tables) ---------------------------- #

def write_count_table(counts: Counter, stats: Dict[str, int], flank: int,
                      include_ambig: bool, collapse_rc: bool, out_path: str) -> None:
    for key in all_context_keys(flank, collapse_rc):
        counts.setdefault(key, 0)
    total = sum(counts.values()) or 1
    headers = ['context', 'left', 'right', 'count', 'fraction']
    with open(out_path, 'w') as out:
        out.write('\t'.join(headers) + '\n')
        def sort_key(k: str):
            left = k[:flank]
            right = k[-flank:] if flank > 0 else ''
            return (left, right, k)
        for ctx in sorted(counts.keys(), key=sort_key):
            left = ctx[:flank]
            right = ctx[-flank:] if flank > 0 else ''
            cnt = counts[ctx]
            frac = cnt / total
            out.write(f"{ctx}\t{left}\t{right}\t{cnt}\t{frac:.6f}\n")
    meta_path = out_path + '.meta.txt'
    with open(meta_path, 'w') as m:
        m.write('CpG Context Counter run stats\n')
        m.write(f"flank\t{flank}\n")
        m.write(f"include_ambiguous\t{include_ambig}\n")
        m.write(f"collapse_revcomp\t{collapse_rc}\n")
        for k in ('total_cpg_sites', 'kept', 'skipped_edge', 'skipped_ambiguous'):
            m.write(f"{k}\t{stats.get(k, 0)}\n")

def append_sites_summary(out_path: str,
                         site_counts: Counter,
                         site_stats: Dict[str,int],
                         flank: int,
                         collapse_rc: bool) -> None:
    for key in all_context_keys(flank, collapse_rc):
        site_counts.setdefault(key, 0)
    total = sum(site_counts.values()) or 1
    headers = ['context', 'left', 'right', 'count', 'fraction']
    with open(out_path, 'a') as out:
        out.write('\n# Per-site contexts (summary from --sites)\n')
        out.write('\t'.join(headers) + '\n')
        def sort_key(k: str):
            left = k[:flank]
            right = k[-flank:] if flank > 0 else ''
            return (left, right, k)
        for ctx in sorted(site_counts.keys(), key=sort_key):
            left = ctx[:flank]
            right = ctx[-flank:] if flank > 0 else ''
            cnt = site_counts[ctx]
            frac = cnt / total
            out.write(f"{ctx}\t{left}\t{right}\t{cnt}\t{frac:.6f}\n")
        out.write('# Per-site summary stats\n')
        for k in ('total_sites','kept_ok','ref_not_CG','out_of_bounds','chrom_missing'):
            out.write(f"# {k}\t{site_stats.get(k, 0)}\n")

# ---------------------------- Main ---------------------------- #

def main():
    p = argparse.ArgumentParser(
        description='Count CpG contexts; optionally add per-site summary, write detailed per-site contexts, and Fisher tests (context + left/right aggregates).'
    )
    p.add_argument('fasta', help='Genome FASTA file (can be multi-FASTA).')
    p.add_argument('-f', '--flank', type=int, default=1, choices=[0, 1, 2, 3],
                   help='Flanking bases on EACH side. 1 => 16 contexts; 2 => 256 contexts.')
    p.add_argument('-o', '--out', default='cpg_contexts.tsv', help='Output TSV path.')
    p.add_argument('-a', '--include-ambig', action='store_true',
                   help='Include contexts with ambiguous bases (N,R,Y, etc.) in the genome-wide counts. Default: skip them.')
    p.add_argument('-r', '--collapse-rc', action='store_true',
                   help='Collapse reverse-complement-equivalent contexts in tables.')
    p.add_argument('--sites', default=None,
                   help='Optional sites table (TSV/CSV/whitespace) with headers: chr, pos, strand. Used for the per-site summary (and detailed file if requested).')
    p.add_argument('--sites-detailed-out', default=None,
                   help='Optional path to write detailed per-site contexts (separate file).')
    p.add_argument('--fishers-out', default='fishers_results.tsv',
                   help='Path to write Fisher test results (default: fishers_results.tsv). Only produced if --sites is provided.')
    args = p.parse_args()

    # Load sites if provided
    site_rows: List[dict] = []
    needed_chroms: Set[str] = set()
    if args.sites:
        site_rows, needed_chroms = load_sites(args.sites)
        if not site_rows:
            print(f"[warn] --sites provided but no rows found in {args.sites}", file=sys.stderr)

    # Genome-wide scan
    grand_counts: Counter = Counter()
    grand_stats = defaultdict(int)
    bg_sides_total = {'left': Counter(), 'right': Counter()}
    chrom_seqs: Dict[str, str] = {}

    for header, seq in parse_fasta(args.fasta):
        counts, stats, side_counts = scan_sequence(seq, flank=args.flank,
                                                   include_ambig=args.include_ambig,
                                                   collapse_rc=args.collapse_rc)
        grand_counts.update(counts)
        for k, v in stats.items():
            grand_stats[k] += v
        for side in ('left','right'):
            bg_sides_total[side].update(side_counts[side])
        if needed_chroms and header in needed_chroms:
            chrom_seqs[header] = seq

    # Write top table
    write_count_table(grand_counts, grand_stats, args.flank,
                      args.include_ambig, args.collapse_rc, args.out)

    # Per-site summary & detailed + Fisher
    if site_rows:
        site_counts, site_stats, fg_sides_total = summarize_sites_contexts(site_rows, chrom_seqs,
                                                                           flank=args.flank,
                                                                           collapse_rc=args.collapse_rc)
        append_sites_summary(args.out, site_counts, site_stats,
                             flank=args.flank, collapse_rc=args.collapse_rc)

        # Fisher tests: per-context + left/right aggregates
        if sum(site_counts.values()) > 0:
            run_fishers(site_counts, grand_counts,
                        flank=args.flank,
                        collapse_rc=args.collapse_rc,
                        fishers_out=args.fishers_out,
                        fg_sides=fg_sides_total,
                        bg_sides=bg_sides_total)
            print(f"Wrote Fisher test results to {args.fishers_out}")
        else:
            print("[warn] No valid per-site contexts to test; Fisher results skipped.", file=sys.stderr)

        if args.sites_detailed_out:
            if not chrom_seqs and needed_chroms:
                print("[warn] No matching chromosomes from --sites found in FASTA; detailed file may be empty/NA.", file=sys.stderr)
            write_sites_detailed(args.sites_detailed_out, site_rows, chrom_seqs, flank=args.flank)

    print(f"Done. Wrote {args.out}" + (f", {args.sites_detailed_out}" if args.sites_detailed_out else ""))
    if site_rows:
        print(f"...and Fisher results: {args.fishers_out}")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
