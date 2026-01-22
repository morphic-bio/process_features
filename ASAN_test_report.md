# ASAN Build Test Report

## Test Setup

Both binaries were built with AddressSanitizer (ASAN) and run with identical parameters on the A375 1k L001 pair dataset.

### Build Commands

**A) PF (khash) ASAN Build:**
```bash
make -C /mnt/pikachu/process_features clean
make CFLAGS="-O1 -g -fsanitize=address -fno-omit-frame-pointer -Iinclude -fopenmp $(pkg-config --cflags cairo)" \
     LDFLAGS="-fsanitize=address -lm -lpthread -lz -fopenmp $(pkg-config --libs cairo) -lhts" \
     assignBarcodes
```

**B) STAR (GLib) ASAN Build:**
```bash
make -C /mnt/pikachu/STAR-suite/core/features/feature_barcodes clean
make CFLAGS="-O1 -g -fsanitize=address -fno-omit-frame-pointer -Iinclude -fopenmp $(pkg-config --cflags glib-2.0) $(pkg-config --cflags cairo)" \
     LDFLAGS="-fsanitize=address -lm -lpthread -lz -fopenmp $(pkg-config --libs glib-2.0) $(pkg-config --libs cairo) -lhts" \
     assignBarcodes
```

### Run Command (both)
```bash
./assignBarcodes \
  --whitelist /storage/A375/3M-5pgex-jan-2023.txt \
  --featurelist /storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv \
  --directory <OUT_DIR> \
  --barcode_fastq_pattern _R1_ --forward_fastq_pattern _R2_ \
  --barcode_length 16 --umi_length 12 --maxHammingDistance 1 \
  --feature_n 1 --barcode_n 2 --max_barcode_mismatches 5 \
  --feature_constant_offset 0 --limit_search -1 --stringency 1 --min_counts 0 \
  --filtered_barcodes /tmp/a375_glib_replacement_test/filtered_barcodes.clean.txt \
  --threads 1 --search_threads 4 --consumer_threads_per_set 1 \
  /storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr/downsampled_1000_L001_pair
```

## Results

### A) PF (khash) ASAN Build
- **Build Status**: ✅ SUCCESS
- **Run Status**: ✅ Completed (exit code 1 due to leaks)
- **Output Files**: ✅ Generated successfully
- **Memory Leaks Detected**: ✅ YES
  - **Summary**: 8376 bytes leaked in 480 allocations
  - **Main leak**: 8075 bytes in 475 objects from `read_barcodes_into_hash` (filtered barcodes)
  - **Other leaks**: File path strings in `organize_fastq_files_by_directory`
- **DEDUP_TRACE Lines**: 0 (feature 11 did not win during dedup for target barcodes)

### B) STAR (GLib) ASAN Build
- **Build Status**: ✅ SUCCESS
- **Run Status**: ✅ Completed (exit code 1 due to leaks)
- **Output Files**: ✅ Generated successfully
- **Memory Leaks Detected**: ✅ YES
  - **Summary**: 301 bytes leaked in 5 allocations
  - **Leaks**: File path strings in `organize_fastq_files_by_directory`
- **DEDUP_TRACE Lines**: 0 (feature 11 did not win during dedup for target barcodes)

### C) Output Comparison

**Filtered Output Statistics:**
- PF barcodes: 475
- STAR barcodes: 475
- Common feature+barcode entries: 479
- Only PF: 1 entry
- Only STAR: 0 entries (but see mismatch below)
- Value differences on common entries: 0

**Specific Mismatch:**
- **PF has**: `Non_Target-1_HD_v2` @ `GTGCAGTTCAGTTAGT` = **479 counts**
- **STAR has**: `Non_Target-1_HD_v2` @ `CAGTAGCTCAATGTTA` = **479 counts**

This is the **same mismatch** as seen in non-ASAN runs: 479 counts of feature 11 (`Non_Target-1_HD_v2`) are assigned to different barcodes.

## Analysis

### Memory Leaks
Both builds have memory leaks, but they are **different**:
- **PF**: Larger leak (8376 bytes) primarily from filtered barcodes hash not being freed
- **STAR**: Smaller leak (301 bytes) from file path strings

The PF leak is in `read_barcodes_into_hash` where `strdup`'d barcode strings are inserted into `khash_t(strptr)` but not freed when the hash is destroyed. This should be fixed by adding proper cleanup.

### DEDUP_TRACE Analysis
**No DEDUP_TRACE lines were logged**, which means:
- Feature 11 (`Non_Target-1_HD_v2`) did **not** win during deduplication for barcodes `CAGTAGCTCAATGTTA` or `GTGCAGTTCAGTTAGT`
- The 479 counts are being assigned **before** deduplication (likely during barcode rescue or initial feature assignment)
- The mismatch occurs in the **barcode rescue** step (`find_best_posterior_match`), not in deduplication

### Root Cause Hypothesis
The mismatch is likely due to **non-deterministic ordering** in barcode rescue when multiple candidates have equal posterior probabilities. The `find_best_posterior_match` function returns the **first** candidate above the threshold, and hash table iteration order differs between `GHashTable` and `khash`.

## Logs Location
- PF ASAN log: `/tmp/a375_pf_asan.log`
- STAR ASAN log: `/tmp/a375_star_asan.log`
- PF output: `/tmp/a375_pf_asan_out/downsampled_1000_L001_pair/filtered/`
- STAR output: `/tmp/a375_star_asan_out/downsampled_1000_L001_pair/filtered/`

## Next Steps
1. **Fix memory leaks**: Add proper cleanup for filtered barcodes hash in PF version
2. **Investigate barcode rescue**: Add logging in `find_best_posterior_match` to capture all candidates and their posteriors for the two problematic barcodes
3. **Make rescue deterministic**: Sort candidates by barcode key when posteriors are equal to ensure deterministic selection
