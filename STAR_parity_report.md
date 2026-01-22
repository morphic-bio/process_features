# STAR Parity Report

## Status: RESOLVED ✅

The GLib-free (khash) build produces **identical output** to the original GLib-based STAR build.

## Summary

After replacing GLib data structures with khash throughout the codebase, comprehensive testing confirmed:

- **475 barcodes** (identical)
- **479 matrix entries** (identical)  
- **787 total counts** (identical)
- **0 value differences** between builds

## Test Configuration

```bash
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
  --filtered_barcodes <filtered_barcodes.txt> \
  --threads 1 \
  --search_threads 4 \
  --consumer_threads_per_set 1 \
  /storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr/downsampled_1000_L001_pair
```

## Changes Made

1. **Replaced GLib with khash**
   - `GHashTable` → `khash_t(...)` hash tables
   - `GArray` → custom `vec_u32_t` and `vec_ptr_t` dynamic arrays
   - `GList` → custom linked list implementations
   - GLib memory functions → standard C `malloc`/`free`

2. **Added khash_wrapper.h**
   - Central header for all khash type definitions and macros

3. **Fixed memory leaks**
   - Added `free_strptr_hash()` for proper cleanup of string-keyed hashes
   - Fixed file path string leaks in `organize_fastq_files_by_directory`
   - Fixed filtered barcodes hash cleanup

## Verification

Both ASAN builds were tested:
- **khash build**: Exit code 0, no application memory leaks
- **GLib build**: Exit code 1, 301 bytes leaked (file paths)

The khash version now has better memory management than the original GLib version.
