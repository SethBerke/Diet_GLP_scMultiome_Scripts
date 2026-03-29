#!/usr/bin/env python3

######################################################################################
# Computes fragment size distributions for all 4 NK samples                          #
# Restricted to GEX-filtered barcodes only                                           #
# Produces nucleosomal banding plot similar to ATAC-Seq library QC figure            #
######################################################################################

import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE        = '/Genomics/pritykinlab/seth/Diet_WL_scMultiome'
CR_BASE     = f'{BASE}/cr_arc_outputs'
BARCODE_DIR = f'{BASE}/Diet_GLP_scMultiome_Scripts/gex_filtered/barcodes'
OUT_DIR     = f'{BASE}/Diet_GLP_scMultiome_Scripts/atac_filtered'

SAMPLES = {
    'NK_SFD': 'NK_1-SFD',
    'NK_HFD': 'NK_2-HFD',
    'NK_GLP': 'NK_3-HFD_GLP',
    'NK_CR':  'NK_4-HFD_CR',
}

MAX_SIZE = 1000  # bp — captures NFR, mono, di, trinucleosome peaks


def load_whitelist(path):
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


def fragment_size_distribution(fragment_path, whitelist, max_size=MAX_SIZE):
    # counts fragment lengths for whitelisted barcodes only
    size_counts = np.zeros(max_size + 1, dtype=np.int64)
    with gzip.open(fragment_path, 'rt') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            fields   = line.rstrip('\n').split('\t')
            barcode  = fields[3]
            if barcode not in whitelist:
                continue
            frag_len = int(fields[2]) - int(fields[1])
            if frag_len <= max_size:
                size_counts[frag_len] += 1
    return size_counts


# collect distributions for all samples
sizes         = np.arange(MAX_SIZE + 1)
distributions = {}

for short_name, cr_name in SAMPLES.items():
    print(f"Processing {short_name}...")
    whitelist     = load_whitelist(f'{BARCODE_DIR}/{short_name}_barcodes.txt')
    fragment_path = f'{CR_BASE}/{cr_name}/outs/atac_fragments.tsv.gz'
    counts        = fragment_size_distribution(fragment_path, whitelist)
    # normalize to relative frequency
    distributions[short_name] = counts / counts.sum()
    print(f"  Done — {counts.sum():,} fragments counted")


def smooth(arr, window=5):
    # rolling average smoothing for cleaner plot
    return pd.Series(arr).rolling(window, center=True).mean().values


# plot
colors = {
    'NK_SFD': 'steelblue',
    'NK_HFD': 'salmon',
    'NK_GLP': 'seagreen',
    'NK_CR':  'mediumpurple'
}

fig, ax = plt.subplots(figsize=(12, 5))

for short_name, dist in distributions.items():
    ax.plot(sizes, smooth(dist), label=short_name,
            color=colors[short_name], linewidth=1.5)

# annotate nucleosomal regions
ax.axvspan(0,   147, alpha=0.05, color='blue',  label='NFR (<147bp)')
ax.axvspan(147, 294, alpha=0.05, color='green', label='Mono-nuc (147-294bp)')
ax.axvspan(294, 500, alpha=0.05, color='red',   label='Di-nuc (294-500bp)')

ax.axvline(147, c='grey', linestyle='--', linewidth=0.8)
ax.axvline(294, c='grey', linestyle='--', linewidth=0.8)
ax.axvline(500, c='grey', linestyle='--', linewidth=0.8)

ax.set_xlabel('Fragment Size (bp)', fontsize=12)
ax.set_ylabel('Relative Frequency', fontsize=12)
ax.set_title('ATAC-Seq Fragment Size Distribution\n(GEX-filtered NK cells)',
             fontsize=13, fontweight='bold')
ax.set_xlim(0, MAX_SIZE)
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/fragment_size_distribution.png', dpi=150, bbox_inches='tight')
plt.show()
print(f"Saved → {OUT_DIR}/fragment_size_distribution.png")