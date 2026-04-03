#!/usr/bin/env python3

import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bisect import bisect_right

BASE        = '/Genomics/pritykinlab/seth/Diet_WL_scMultiome'
CR_BASE     = f'{BASE}/cr_arc_outputs'
BARCODE_DIR = f'{BASE}/Diet_GLP_scMultiome_Scripts/gex_filtered/barcodes'
OUT_DIR     = f'{BASE}/Diet_GLP_scMultiome_Scripts/atac_filtered'
TSS_BED     = f'{BASE}/mm10_tss.bed'

SAMPLES = {
    'NK_SFD': 'NK_1-SFD',
    'NK_HFD': 'NK_2-HFD',
    'NK_GLP': 'NK_3-HFD_GLP',
    'NK_CR':  'NK_4-HFD_CR',
}

WINDOW = 2000

def load_whitelist(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())

def load_tss(tss_bed):
    tss = pd.read_csv(tss_bed, sep='\t', header=None,
                      names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    tss = tss.drop_duplicates(subset=['chr', 'start'])
    tss_lookup = {}
    for chrom, sub in tss.groupby('chr', sort=False):
        tss_lookup[chrom] = np.sort(sub['start'].to_numpy(dtype=np.int64))
    return tss_lookup

def compute_aggregate_tss_profile(fragment_path, whitelist, tss_lookup, window=2000):
    n_bins  = 2 * window + 1
    profile = np.zeros(n_bins, dtype=np.float64)
    with gzip.open(fragment_path, 'rt') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            fields  = line.rstrip('\n').split('\t')
            chrom, start_text, end_text, barcode = fields[:4]
            if barcode not in whitelist:
                continue
            chrom_tss = tss_lookup.get(chrom)
            if chrom_tss is None:
                continue
            for pos in (int(start_text) + 4, int(end_text) - 5):
                idx_left  = bisect_right(chrom_tss, pos + window)
                idx_right = bisect_right(chrom_tss, pos - window)
                for tss_pos in chrom_tss[idx_right:idx_left]:
                    offset = pos - tss_pos + window
                    if 0 <= offset < n_bins:
                        profile[offset] += 1
    flank      = 100
    background = (profile[:flank].mean() + profile[-flank:].mean()) / 2
    profile    = profile / max(background, 1e-6)
    return profile

print("Loading TSS...")
tss_lookup = load_tss(TSS_BED)

positions = np.arange(-WINDOW, WINDOW + 1)
profiles  = {}
colors    = {'NK_SFD': 'steelblue', 'NK_HFD': 'salmon',
             'NK_GLP': 'seagreen',  'NK_CR':  'mediumpurple'}

for short_name, cr_name in SAMPLES.items():
    print(f"Processing {short_name}...")
    whitelist     = load_whitelist(f'{BARCODE_DIR}/{short_name}_barcodes.txt')
    fragment_path = f'{CR_BASE}/{cr_name}/outs/atac_fragments.tsv.gz'
    profiles[short_name] = compute_aggregate_tss_profile(
        fragment_path, whitelist, tss_lookup, WINDOW)
    print(f"  Done — peak enrichment: {profiles[short_name][WINDOW]:.2f}x")

def smooth(arr, window=10):
    return pd.Series(arr).rolling(window, center=True).mean().values

fig, ax = plt.subplots(figsize=(10, 6))
for short_name, profile in profiles.items():
    ax.plot(positions, smooth(profile), label=short_name,
            color=colors[short_name], linewidth=1.5)

ax.axvline(0,    c='grey', linestyle='--', linewidth=0.8)
ax.axhline(1,    c='grey', linestyle=':',  linewidth=0.8)
ax.axvspan(-147, 147, alpha=0.05, color='blue', label='NFR')
ax.set_xlabel('Distance From TSS (bp)', fontsize=12)
ax.set_ylabel('Normalized Insertion Profile', fontsize=12)
ax.set_title('TSS Enrichment Profile\n(GEX-filtered NK cells)',
             fontsize=13, fontweight='bold')
ax.set_xlim(-WINDOW, WINDOW)
ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/tss_enrichment_profile.png', dpi=150, bbox_inches='tight')
print(f"Saved → {OUT_DIR}/tss_enrichment_profile.png")