# STAR Parity Investigation Report

## Context
User requirement: **New GLib-free build must be identical to old STAR results** (Cell Ranger parity can be slightly different, but STAR-to-STAR must match exactly).

I compared a baseline STAR output against the new GLib-free output and found a single barcode-level discrepancy even though totals match.

## Baseline vs New Outputs

**Baseline (old STAR build)**
- Directory: `/tmp/a375_cb_rescue_scan_1000_L001_pair_m1_mbm5/m1_fn1_bn2_mbm5_off0_lsn1/downsampled_1000_L001_pair/filtered`
- Counts summary:
  - `barcodes.txt`: 475
  - `matrix.mtx` total counts: 787

**New (GLib-free build)**
- Directory: `/tmp/a375_glib_replacement_test_clean3/downsampled_1000_L001_pair/filtered`
- Counts summary:
  - `barcodes.txt`: 475
  - `matrix.mtx` total counts: 787

Totals and barcode sets match.

## Specific Mismatch
Only one **feature+barcode** entry differs:

- Baseline has:
  - `Non_Target-1_HD_v2` @ `CAGTAGCTCAATGTTA` = **479**
- New has:
  - `Non_Target-1_HD_v2` @ `GTGCAGTTCAGTTAGT` = **479**

Everything else matches exactly. This is a **swap of 479 counts** between two barcodes for the same feature.

Per-feature/per-barcode snapshot:
- Baseline:
  - `CAGTAGCTCAATGTTA`: `Non_Target-1_HD_v2=479`, `RAB1A-2_MS=3`, `Non_Target-1_MS=2`
  - `GTGCAGTTCAGTTAGT`: `Non_Target-1_MS=1`
- New:
  - `CAGTAGCTCAATGTTA`: `RAB1A-2_MS=3`, `Non_Target-1_MS=2`
  - `GTGCAGTTCAGTTAGT`: `Non_Target-1_HD_v2=479`, `Non_Target-1_MS=1`

So the **479 counts of `Non_Target-1_HD_v2` moved between those two barcodes**.

## Likely Area
This looks like **barcode rescue tie-breaking/order** rather than feature assignment:
- Feature list is identical.
- Barcode set is identical.
- Only one barcode pair is affected.

Likely in `assignBarcodes.c` around barcode rescue logic:
- `find_closest_barcodes`
- `find_best_posterior_match`
- `process_pending_barcodes`

Specifically, iteration/ordering differences between `GHashTable` and `khash` could yield a different rescue decision when multiple closest candidates are possible.

## Reproduction
Command used (both baseline and new runs):
```
./assignBarcodes \
  --whitelist /storage/A375/3M-5pgex-jan-2023.txt \
  --featurelist /storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv \
  --directory <OUT_DIR> \
  --barcode_fastq_pattern _R1_ \
  --forward_fastq_pattern _R2_ \
  --barcode_length 16 \
  --umi_length 12 \
  --maxHammingDistance 1 \
  --feature_n 1 \
  --barcode_n 2 \
  --max_barcode_mismatches 5 \
  --feature_constant_offset 0 \
  --limit_search -1 \
  --stringency 1 \
  --min_counts 0 \
  --filtered_barcodes /tmp/a375_glib_replacement_test/filtered_barcodes.clean.txt \
  --threads 1 \
  --search_threads 4 \
  --consumer_threads_per_set 1 \
  /storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr/downsampled_1000_L001_pair
```

`filtered_barcodes.clean.txt` is the CR barcode list with `-1` suffix stripped:
```
sed 's/-1$//' /tmp/a375_crispr_1000_pair/outs/multi/count/raw_feature_bc_matrix_crispr_only_nonzero/barcodes.tsv > /tmp/a375_glib_replacement_test/filtered_barcodes.clean.txt
```

## Next Steps for Agent
1. Add debug logging in barcode rescue:
   - Log candidates + posteriors for any rescue involving:
     - `CAGTAGCTCAATGTTA`
     - `GTGCAGTTCAGTTAGT`
2. Compare against baseline if possible (or re-run baseline build with equivalent logging).
3. Identify deterministic tie-breaking differences in:
   - `find_best_posterior_match`
   - `process_pending_barcodes`
   - any iteration across khash entries affecting ordering
4. Fix to ensure **byte-for-byte identical** rescue selection to STAR baseline.

