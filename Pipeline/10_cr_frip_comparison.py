#!/usr/bin/env python3

######################################################################################
# Computes insertion-based FRiP against Cell Ranger called peaks                     #
# Uses same binary search method as atlas FRiP for apples-to-apples comparison       #
# Outputs per-sample median FRiP scores only                                         #
######################################################################################

import gzip
import os
import numpy as np
import pandas as pd
from bisect import bisect_right

BASE        = '/Genomics/pritykinlab/seth/Diet_WL_scMultiome'
CR_BASE     = f'{BASE}/cr_arc_outputs'
BARCODE_DIR = f'{BASE}/Diet_GLP_scMultiome_Scripts/gex_filtered/barcodes'
ATAC_DIR    = f'{BASE}/Diet_GLP_scMultiome_Scripts/atac_filtered'

SAMPLES = {
    'NK_SFD': 'NK_1-SFD',
    'NK_HFD': 'NK_2-HFD',
    'NK_GLP': 'NK_3-HFD_GLP',
    'NK_CR':  'NK_4-HFD_CR',
}


def load_whitelist(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


def load_cr_peaks(peaks_bed):
    # reads CR peaks BED (3 columns: chr start end) and builds per-chromosome sorted arrays
    peaks = pd.read_csv(peaks_bed, sep="\t", header=None,
                        names=["chr", "start", "end"],
                        comment="#")
    chrom_lookup = {}
    for chrom, sub in peaks.groupby("chr", sort=False):
        sub    = sub.sort_values("start")
        starts = sub["start"].to_numpy(dtype=np.int64)
        ends   = sub["end"].to_numpy(dtype=np.int64)
        chrom_lookup[chrom] = (starts, ends)
    return chrom_lookup, len(peaks)


def compute_cr_frip(fragment_path, whitelist, chrom_lookup):
    # streams fragment file and computes insertion-based FRiP against CR peaks
    # numerator:   Tn5 insertion events (start+4, end-5) overlapping CR called peaks
    # denominator: total Tn5 insertion events across genome (2 per fragment)
    barcode_to_idx  = {b: i for i, b in enumerate(sorted(whitelist))}
    n_cells         = len(barcode_to_idx)
    total_ins       = np.zeros(n_cells, dtype=np.int64)
    peak_ins        = np.zeros(n_cells, dtype=np.int64)

    with gzip.open(fragment_path, "rt") as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields  = line.rstrip("\n").split("\t")
            chrom, start_text, end_text, barcode = fields[:4]
            cell_idx = barcode_to_idx.get(barcode)
            if cell_idx is None:
                continue

            total_ins[cell_idx] += 2

            chrom_data = chrom_lookup.get(chrom)
            if chrom_data is None:
                continue

            starts, ends = chrom_data
            for pos in (int(start_text) + 4, int(end_text) - 5):
                candidate = bisect_right(starts, pos) - 1
                if candidate >= 0 and pos < ends[candidate]:
                    peak_ins[cell_idx] += 1

    frip = peak_ins / np.maximum(total_ins, 1)
    return frip


# compute and compare FRiP scores
records = []
for short_name, cr_name in SAMPLES.items():
    print(f"Processing {short_name}...")

    whitelist     = load_whitelist(f'{BARCODE_DIR}/{short_name}_barcodes.txt')
    fragment_path = f'{CR_BASE}/{cr_name}/outs/atac_fragments.tsv.gz'
    cr_peaks_bed  = f'{CR_BASE}/{cr_name}/outs/atac_peaks.bed'

    chrom_lookup, n_peaks = load_cr_peaks(cr_peaks_bed)
    frip_scores           = compute_cr_frip(fragment_path, whitelist, chrom_lookup)

    records.append({
        'Sample':        short_name,
        'N Cells':       len(whitelist),
        'CR Peaks':      n_peaks,
        'CR FRiP (ins)': round(np.median(frip_scores), 3),
    })
    print(f"  CR peaks: {n_peaks:,} | Median FRiP: {np.median(frip_scores):.3f}")

# print final comparison table
print("\n=== FRiP Comparison ===")
df = pd.DataFrame(records)
print(df.to_string(index=False))