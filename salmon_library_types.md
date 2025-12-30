# Salmon Library Type: ISR vs Auto Mode

## Overview

Salmon and related tools (alevin, alevin-fry, piscem) use library type specifications to understand the orientation and strandedness of sequencing reads. This is critical for accurate quantification in RNA-seq experiments.

## Library Type Format

In salmon, library types follow the format: `<STRANDEDNESS><ORIENTATION><DIRECTION>`

For single-cell RNA-seq with alevin/simpleaf, the commonly used parameter is `--expected-ori` which simplifies this specification.

## ISR Mode (Inward Stranded Reverse)

**ISR** stands for **Inward, Stranded, Reverse-complement**.

### Characteristics:
- **I (Inward)**: Paired-end reads face inward toward each other (most common for Illumina paired-end sequencing)
- **S (Stranded)**: The library is strand-specific (preserves information about which DNA strand was transcribed)
- **R (Reverse)**: The first read (R1) comes from the reverse complement of the transcript

### Behavior:
- Salmon explicitly expects reads to match the ISR orientation pattern
- Reads that don't match this pattern may be incorrectly quantified or discarded
- This is the standard for many stranded RNA-seq protocols including:
  - Illumina TruSeq Stranded mRNA
  - 10x Chromium single-cell RNA-seq (uses `fw` which is equivalent for single-end barcode reads)
  - dUTP-based stranded protocols

### Use Case:
Use ISR when you **know** your library is stranded and follows the reverse-complement orientation. This provides:
- ✅ More accurate quantification
- ✅ Strand-specific information
- ✅ Better handling of overlapping genes on opposite strands
- ❌ Incorrect results if library prep doesn't match

## Auto Mode (Automatic Detection)

**Auto** (or **A**) enables automatic library type detection.

### Characteristics:
- Salmon examines a subset of reads at the beginning of quantification
- Infers the most likely library type from the actual data
- Falls back to the most probable orientation if uncertain

### Behavior:
- Analyzes the first ~5 million reads (configurable)
- Tests different orientation patterns to find the best match
- Reports the detected library type in the output
- Can handle mixed or uncertain library preparations

### Use Case:
Use auto mode when:
- ✅ Library type is unknown or uncertain
- ✅ Testing new protocols or library prep methods
- ✅ Working with data from external sources
- ✅ Prefer robustness over strict assumptions
- ❌ Slightly slower than explicit mode (due to detection phase)
- ⚠️ May incorrectly detect type with low-quality or unusual data

## Key Differences Summary

| Aspect | ISR Mode | Auto Mode |
|--------|----------|-----------|
| **Speed** | Faster (no detection phase) | Slightly slower (detection overhead) |
| **Accuracy** | Highest when library matches | Good, adapts to actual data |
| **Robustness** | Fails if library doesn't match | Handles uncertainty |
| **Use Case** | Known, standard protocols | Unknown or variable protocols |
| **Error Handling** | Silent incorrect quantification | May detect unexpected patterns |

## Current Usage in This Repository

In `scripts/runMultiAlign.sh`, the repository uses:

```bash
--expected-ori fw
```

For single-cell RNA-seq with 10x Chromium:
- `fw` (forward) is appropriate because the biological read (R2) is in the forward orientation
- R1 contains the cell barcode + UMI (technical read)
- R2 contains the actual cDNA sequence (biological read)

This is functionally similar to specifying the strand orientation for the quantification step.

## Recommendations

### Use ISR (or fw for single-cell) when:
1. Working with **standard 10x Chromium single-cell RNA-seq**
2. Library prep protocol is **well-documented** and confirmed
3. Processing **large production datasets** where speed matters
4. Protocol has been **validated** on test samples

### Use Auto mode when:
1. **Protocol is uncertain** or from external sources
2. **Testing new library prep methods**
3. **Troubleshooting** quantification issues
4. Working with **mixed samples** or unusual prep
5. Want to **verify** the library orientation assumption

## Checking Library Type

After running salmon/alevin with auto mode, check the log file for:
```
Library format: ISR
```

This confirms what salmon detected and can inform future runs with explicit mode.

## References

- [Salmon Documentation - Library Type](https://salmon.readthedocs.io/en/latest/library_type.html)
- [Alevin-fry Documentation](https://alevin-fry.readthedocs.io/)
- [10x Genomics Library Structure](https://www.10xgenomics.com/)

## Impact on This Pipeline

For the `runMultiAlign.sh` script:
- Current implementation uses `--expected-ori fw` (explicit mode)
- This is correct for standard 10x single-cell RNA-seq
- To enable auto-detection, change to `--expected-ori both` or omit the parameter (depending on tool version)

**Example modification for auto mode:**
```bash
# Current (explicit):
--expected-ori fw --unfiltered-pl

# Auto mode option 1:
--expected-ori both --unfiltered-pl

# Auto mode option 2 (if supported):
--expected-ori auto --unfiltered-pl
```

Note: The specific syntax may vary between salmon, piscem, and simpleaf versions. Always consult the tool's documentation for the version you're using.
