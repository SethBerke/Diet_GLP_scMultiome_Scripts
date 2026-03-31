#!/usr/bin/env python3

######################################################################################
# Builds ATAC AnnData object for a single sample                                     #
# Scored against the NK Atlas Peak Set using binary search on Tn5 insertion sites    #
# Computes FRiP, TSS enrichment, and nucleosomal signal per cell                     #
######################################################################################

import argparse
import gzip
import os
import subprocess
import tempfile
from bisect import bisect_right

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ATLAS_COLUMNS = ["chr", "start", "end", "peak_id", "source", "type"]


def condition_from_sample(sample):
    # extracts condition from sample name e.g. NK_SFD -> SFD
    return sample.split("_", 1)[1]

# reads atlas BED file and builds per-chromosome sorted arrays for binary search
def load_atlas(atlas_bed):
    atlas = pd.read_csv(atlas_bed, sep="\t", header=None, names=ATLAS_COLUMNS)
    atlas["interval"] = atlas["chr"] + ":" + atlas["start"].astype(str) + "-" + atlas["end"].astype(str)
    chrom_lookup = {}
    for chrom, sub in atlas.groupby("chr", sort=False):
        sub     = sub.sort_values("start")  # ensure sorted within chromosome
        starts  = sub["start"].to_numpy(dtype=np.int64, copy=True)
        ends    = sub["end"].to_numpy(dtype=np.int64, copy=True)
        indices = sub.index.to_numpy(dtype=np.int64, copy=True)
        chrom_lookup[chrom] = (starts, ends, indices)
    return atlas, chrom_lookup


def load_whitelist(whitelist_path):
    # reads GEX-filtered barcodes — only these cells will be processed
    with open(whitelist_path) as f:
        return set(line.strip() for line in f if line.strip())


def load_tss(tss_bed):
    # reads TSS BED file and builds per-chromosome sorted arrays of TSS positions
    tss = pd.read_csv(tss_bed, sep="\t", header=None,
                      names=["chr", "start", "end", "name", "score", "strand"])
    tss = tss.drop_duplicates(subset=["chr", "start"])
    tss_lookup = {}
    for chrom, sub in tss.groupby("chr", sort=False):
        tss_lookup[chrom] = np.sort(sub["start"].to_numpy(dtype=np.int64))
    return tss_lookup


def write_mapped_insertions(fragment_path, barcode_to_idx, chrom_lookup, mapped_path):
    # streams fragment file once; for each whitelisted barcode:
    #   - counts total insertions (2 per fragment, weight=1 always)
    #   - counts sub-nucleosomal (<147bp) and mono-nucleosomal (147-294bp) fragments
    #   - maps each Tn5 insertion site (start+4, end-5) to atlas peaks via binary search
    #   - writes (cell_idx, peak_idx, 1) to mapped_path for sparse matrix construction
    # note: column 5 of fragment file is duplicate count across barcodes, not a per-cell weight
    n_cells          = len(barcode_to_idx)
    total_insertions = np.zeros(n_cells, dtype=np.int64)
    peak_insertions  = np.zeros(n_cells, dtype=np.int64)
    sub_nuc_counts   = np.zeros(n_cells, dtype=np.int64)
    mono_nuc_counts  = np.zeros(n_cells, dtype=np.int64)

    with gzip.open(fragment_path, "rt") as handle, open(mapped_path, "w") as out_handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start_text, end_text, barcode = fields[:4]
            cell_idx = barcode_to_idx.get(barcode)
            if cell_idx is None:
                continue

            # always weight=1; column 5 is experiment-wide duplicate count, not per-cell
            total_insertions[cell_idx] += 2

            # fragment length for nucleosomal signal
            frag_len = int(end_text) - int(start_text)
            if frag_len < 147:
                sub_nuc_counts[cell_idx]  += 1
            elif frag_len < 294:
                mono_nuc_counts[cell_idx] += 1

            # map Tn5 insertion sites to atlas peaks
            chrom_data = chrom_lookup.get(chrom)
            if chrom_data is None:
                continue

            starts, ends, indices = chrom_data
            for pos in (int(start_text) + 4, int(end_text) - 5):
                candidate = bisect_right(starts, pos) - 1
                if candidate >= 0 and pos < ends[candidate]:
                    out_handle.write(f"{cell_idx}\t{indices[candidate]}\t1\n")
                    peak_insertions[cell_idx] += 1

    return total_insertions, peak_insertions, sub_nuc_counts, mono_nuc_counts


def sort_mapped_pairs(mapped_path, sorted_path, tmpdir, threads):
    # sorts (cell_idx, peak_idx) pairs by cell then peak for efficient sparse matrix construction
    with open(sorted_path, "w") as out_handle:
        subprocess.run(
            ["sort", "--temporary-directory", tmpdir,
             "--parallel", str(threads), "-k1,1n", "-k2,2n", mapped_path],
            check=True, stdout=out_handle,
            env={**os.environ, "LC_ALL": "C"},
        )


def sorted_pairs_to_matrix(sorted_path, n_cells, n_peaks):
    # reads sorted (cell_idx, peak_idx, weight) pairs and accumulates into a sparse CSR matrix
    rows, cols, data = [], [], []
    last_pair, pair_count = None, 0

    with open(sorted_path) as handle:
        for line in handle:
            cell_idx_text, peak_idx_text, weight_text = line.rstrip("\n").split("\t")
            pair   = (int(cell_idx_text), int(peak_idx_text))
            weight = int(weight_text)
            if pair == last_pair:
                pair_count += weight
                continue
            if last_pair is not None:
                rows.append(last_pair[0]); cols.append(last_pair[1]); data.append(pair_count)
            last_pair, pair_count = pair, weight

    if last_pair is not None:
        rows.append(last_pair[0]); cols.append(last_pair[1]); data.append(pair_count)

    return sparse.csr_matrix(
        (np.asarray(data, dtype=np.int32), (np.asarray(rows), np.asarray(cols))),
        shape=(n_cells, n_peaks), dtype=np.int32,
    )


def compute_tsse(fragment_path, barcode_to_idx, tss_lookup, window=2000):
    # computes per-cell TSS enrichment score by building insertion profiles around all TSSs
    # score = mean insertion rate at +/-50bp around TSS / mean insertion rate at outer 100bp flanks
    # higher score = more signal at open chromatin promoter regions relative to background
    n_cells  = len(barcode_to_idx)
    n_bins   = 2 * window + 1
    profiles = np.zeros((n_cells, n_bins), dtype=np.float32)

    with gzip.open(fragment_path, "rt") as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start_text, end_text, barcode = fields[:4]
            cell_idx = barcode_to_idx.get(barcode)
            if cell_idx is None:
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
                        profiles[cell_idx, offset] += 1

    flank      = 100
    background = (profiles[:, :flank].mean(axis=1) +
                  profiles[:, -flank:].mean(axis=1)) / 2
    center     = profiles[:, window - 50: window + 50].mean(axis=1)
    tsse       = center / np.maximum(background, 1e-6)

    return tsse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fragment-file', required=True)
    parser.add_argument('--whitelist',     required=True)
    parser.add_argument('--atlas',         required=True)
    parser.add_argument('--tss-bed',       required=True)
    parser.add_argument('--output',        required=True)
    parser.add_argument('--sample-name',   required=True)
    parser.add_argument('--threads',       type=int, default=4)
    args = parser.parse_args()

    print(f"\n=== Processing {args.sample_name} ===")

    # load all reference data and whitelist
    print("Loading atlas, TSS, and whitelist...")
    atlas, chrom_lookup = load_atlas(args.atlas)
    tss_lookup          = load_tss(args.tss_bed)
    whitelist           = load_whitelist(args.whitelist)
    barcode_to_idx      = {b: i for i, b in enumerate(sorted(whitelist))}
    barcodes            = sorted(whitelist)
    n_cells             = len(barcodes)
    n_peaks             = len(atlas)
    print(f"  Atlas peaks:        {n_peaks:,}")
    print(f"  Whitelisted cells:  {n_cells:,}")
    print(f"  TSS sites:          {sum(len(v) for v in tss_lookup.values()):,}")

    with tempfile.TemporaryDirectory(prefix="atac_build_") as tmpdir:

        mapped_path = os.path.join(tmpdir, "mapped.tsv")
        sorted_path = os.path.join(tmpdir, "sorted.tsv")

        # step 1: stream fragment file → insertion counts + nucleosomal signal + peak matrix
        print("Mapping Tn5 insertions to atlas peaks...")
        total_ins, peak_ins, sub_nuc, mono_nuc = write_mapped_insertions(
            args.fragment_file, barcode_to_idx, chrom_lookup, mapped_path
        )

        # step 2: sort mapped pairs and build sparse cell x peak matrix
        print("Sorting and building peak matrix...")
        sort_mapped_pairs(mapped_path, sorted_path, tmpdir, args.threads)
        matrix = sorted_pairs_to_matrix(sorted_path, n_cells, n_peaks)
        print(f"  Matrix shape: {matrix.shape}")

        # step 3: stream fragment file again to compute per-cell TSS enrichment
        print("Computing TSS enrichment...")
        tsse_scores = compute_tsse(args.fragment_file, barcode_to_idx, tss_lookup)
        print(f"  TSS range: {tsse_scores.min():.2f} - {tsse_scores.max():.2f}")

    # step 4: assemble obs dataframe with all QC metrics
    obs = pd.DataFrame(index=[f"{args.sample_name}:{b}" for b in barcodes])
    obs["barcode"]            = barcodes
    obs["sample"]             = args.sample_name
    obs["condition"]          = condition_from_sample(args.sample_name)
    obs["total_insertions"]   = total_ins
    obs["peak_insertions"]    = peak_ins
    obs["FRiP"]               = peak_ins  / np.maximum(total_ins, 1)
    obs["nucleosomal_signal"] = mono_nuc  / np.maximum(sub_nuc,   1)
    obs["n_counts"]           = matrix.sum(axis=1).A1
    obs["n_features_per_cell"] = np.asarray((matrix > 0).sum(axis=1)).flatten()
    obs["tsse"]               = tsse_scores

    print(f"\n  Median FRiP:               {obs['FRiP'].median():.3f}")
    print(f"  Median TSS:                {obs['tsse'].median():.2f}")
    print(f"  Median nucleosomal_signal: {obs['nucleosomal_signal'].median():.2f}")
    print(f"  Median total_insertions:   {obs['total_insertions'].median():.0f}")
    print(f"  Median n_features_per_cell: {obs['n_features_per_cell'].median():.0f}")


    # step 5: build var dataframe from atlas peaks
    var       = atlas.copy()
    var.index = var["peak_id"].astype(str)

    # step 6: save as h5ad
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.write_h5ad(args.output)
    print(f"\nSaved {args.sample_name}: {adata.shape} -> {args.output}")


if __name__ == "__main__":
    main()