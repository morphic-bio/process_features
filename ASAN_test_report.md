# ASAN Build Test Report

## Status: PASSED ✅

Both binaries were built with AddressSanitizer (ASAN) and tested. The khash version has **no application memory leaks** after fixes.

## Build Commands

**khash Build (GLib-free):**
```bash
make CFLAGS="-O1 -g -fsanitize=address -fno-omit-frame-pointer -Iinclude -fopenmp $(pkg-config --cflags cairo)" \
     LDFLAGS="-fsanitize=address -lm -lpthread -lz -fopenmp $(pkg-config --libs cairo) -lhts" \
     assignBarcodes
```

**GLib Build (original):**
```bash
make CFLAGS="-O1 -g -fsanitize=address -fno-omit-frame-pointer -Iinclude -fopenmp $(pkg-config --cflags glib-2.0) $(pkg-config --cflags cairo)" \
     LDFLAGS="-fsanitize=address -lm -lpthread -lz -fopenmp $(pkg-config --libs glib-2.0) $(pkg-config --libs cairo) -lhts" \
     assignBarcodes
```

## Results

### khash Build (after fixes)
- **Build Status**: ✅ SUCCESS
- **Run Status**: ✅ Exit code 0
- **Memory Leaks**: ✅ None (only ~1KB from external Cairo/fontconfig libraries)

### GLib Build (original)
- **Build Status**: ✅ SUCCESS  
- **Run Status**: ⚠️ Exit code 1 (due to leaks)
- **Memory Leaks**: 301 bytes in 5 allocations (file path strings)

## Output Comparison

| Metric | khash | GLib |
|--------|-------|------|
| Barcodes | 475 | 475 |
| Entries | 479 | 479 |
| Total counts | 787 | 787 |
| Common entries | 479 | 479 |
| Differences | 0 | 0 |

**Outputs are byte-for-byte identical.**

## Memory Leak Fixes Applied

1. **`free_strptr_hash()`** - New function to properly free strdup'd keys in string-keyed hash tables
2. **`free_fastq_files_collection()`** - Updated to free individual file path strings and sample names
3. **`organize_fastq_files_by_directory()`** - Added cleanup of temporary file path arrays

## Conclusion

The khash version:
- Produces identical output to the GLib version
- Has better memory management (no application leaks)
- Removes the GLib dependency for cleaner builds
