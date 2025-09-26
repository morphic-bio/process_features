# Code File Details - Technical Documentation

## Source Files (src/)

### Core Application Files

#### main.c
**Purpose**: Program entry point and command-line interface
**Key Functions**:
- `main()`: Argument parsing, initialization, and process orchestration
- `print_usage()`: Comprehensive CLI help display
- Multi-process fork/wait pattern for parallel sample processing

**Technical Details**:
- Uses `getopt_long()` for extensive CLI parsing (42+ options)
- Initializes global lookup tables: `seq2code`, `code2seq`, `diff2Hamming`
- Creates shared memory for thread counting using `mmap()`
- Forks separate processes for each sample with configurable concurrency limits
- Memory management through `sample_args` structure passed to child processes

**Key Data Structures**:
- `sample_args`: Comprehensive parameter structure for sample processing
- `fastq_files_collection`: File organization and sample mapping
- Global hash tables for whitelist and feature codes

#### assignBarcodes.c
**Purpose**: Core barcode assignment and UMI deduplication logic
**Key Functions**:
- `process_files_in_sample()`: Main sample processing orchestrator
- `find_deduped_counts()`: UMI clique analysis and deduplication
- `add_deduped_count()`: Stringency-based count assignment
- `printFeatureCounts()`: Matrix Market output generation
- `checkAndCorrectBarcode()`: Barcode error correction with posterior probabilities

**Technical Details**:
- Producer-consumer threading model with ring buffer I/O
- Connected component analysis for UMI error correction (Hamming distance ≤1)
- Three stringency modes: RNA-seq (0), proportional (1-999), strict (≥1000)
- Memory pool allocation for all dynamic structures
- Interactive histogram and heatmap generation

**Key Algorithms**:
- Graph-based UMI deduplication using breadth-first search
- Bayesian posterior probability for barcode rescue
- Fuzzy matching with configurable Hamming distance thresholds
- Exhaustive vs. targeted feature search strategies

#### barcode_match.c
**Purpose**: Sequence encoding, feature lookup, and matching utilities
**Key Functions**:
- `string2code()`: Convert DNA sequences to 2-bit packed representation
- `code2string()`: Decode packed sequences back to DNA strings
- `feature_lookup_kmer()`: Optimized k-mer lookup with direct/hash modes
- `check_sequence()`: Validate DNA sequences (ACGT only)

**Technical Details**:
- 2-bit DNA encoding (A=0, C=1, G=2, T=3) for memory efficiency
- Direct 64-bit comparison for small feature sets (≤128 variants)
- Hash table fallback using GLib GBytes for larger feature sets
- Precomputed lookup tables for sequence conversion and Hamming distances

**Optimization Features**:
- SIMD-friendly bit operations for sequence comparison
- Cache-efficient direct array lookup for small probe sets
- Lazy initialization of 64-bit probe tables

#### demux_bam.c
**Purpose**: BAM file processing for STAR Solo output
**Key Functions**:
- `process_bam_single()`: Main BAM processing loop
- `extract_kmer_from_bam()`: Extract probe k-mers from BAM sequence data
- Probe detection with configurable offset and nearby search

**Technical Details**:
- HTSlib integration for BAM reading and tag extraction
- Extracts CB (cell barcode), UB (UMI), GX/GE (gene) tags
- Handles both forward and reverse-strand reads with complement calculation
- Memory-efficient counting without heavy deduplication tables
- Configurable filtering: MAPQ, duplicates, alignment type

**BAM Processing Features**:
- Primary alignment filtering with optional secondary/supplementary inclusion
- Duplicate read filtering via BAM_FDUP flag
- Intergenic read handling (GX='-') with optional inclusion
- Probe search with fallback offset testing (±1, ±2)

#### demux_fastq.c
**Purpose**: FASTQ-based sample demultiplexing
**Key Functions**:
- `process_fastq_single()`: Main FASTQ processing with producer-consumer model
- `extract_probe_and_demux()`: Probe extraction and sample assignment
- `get_or_open_sink()`: Lazy output file creation and management

**Technical Details**:
- Multi-threaded producer-consumer with ring buffer
- Per-sample output sinks with mutex-protected concurrent writing
- Direct 64-bit probe comparison vs. hash table lookup
- Automatic file organization with gzip compression

**Performance Optimizations**:
- Lock-free ring buffer for high-throughput I/O
- Per-sink mutexes to prevent write interleaving
- Configurable buffer sizes and thread counts

#### io.c
**Purpose**: File I/O operations and data structure initialization
**Key Functions**:
- `read_features_file()`: Parse CSV feature files with flexible column detection
- `organize_fastq_files_by_directory()`: Automatic FASTQ file discovery and grouping
- `allocate_feature_arrays()`: Memory allocation for feature data structures
- `read_whiteList()`: Barcode whitelist loading with validation

**Technical Details**:
- Flexible CSV parsing with header detection
- Pattern-based FASTQ file matching (_R1_, _R2_, _R3_)
- Memory-efficient feature array allocation with aligned storage
- Whitelist validation with sequence length checking and error correction

**File Organization Features**:
- Automatic sample detection from directory structure
- Support for explicit file lists vs. directory scanning
- Pattern matching for different sequencing platforms
- Gzip-aware file handling

#### memory.c
**Purpose**: Custom memory pool management system
**Key Functions**:
- `initialize_memory_pool_collection()`: Create pools for all data types
- `allocate_memory_from_pool()`: Fast pool-based allocation
- `expand_memory_pool()`: Dynamic pool expansion
- `free_memory_pool_collection()`: Cleanup all pools

**Technical Details**:
- Pre-allocated memory blocks for different structure types
- Block-based allocation with automatic expansion
- Zero-copy pool expansion with linked block lists
- Thread-safe allocation within individual pools

**Pool Types**:
- `feature_counts_pool`: Barcode-feature count matrices
- `feature_umi_counts_pool`: UMI-level count tracking
- `feature_sequences_pool`: Matched sequence storage
- `unmatched_barcodes_features_block_pool`: Error correction candidates
- `cb_counts_pool`: Cell barcode count arrays

#### heatmap.c
**Purpose**: QC heatmap generation using Cairo graphics
**Key Functions**:
- `generate_heatmap()`: Feature richness heatmap creation
- `generate_deduped_heatmap()`: UMI count distribution heatmap
- `draw_heatmap_core()`: Shared rendering engine
- `select_colormap()`: Dynamic colormap selection based on data range

**Technical Details**:
- Cairo-based PNG generation with anti-aliasing
- Dynamic colormap selection (16/64/256/1024 colors) based on data range
- Integrated bar graphs showing column/row summaries
- Configurable cell sizes and padding for different data scales

**Visual Features**:
- Plasma colormap for perceptually uniform color scaling
- Feature name labels with automatic text sizing
- Color bars with value annotations
- Bar graphs showing marginal distributions

#### plot_histogram.c
**Purpose**: Interactive HTML histogram generation
**Key Functions**:
- `plot_simple_histogram()`: Generate Plotly.js-based interactive plots
- `generate_histogram_html()`: HTML template with embedded data
- Support for linear/log scale switching

**Technical Details**:
- Plotly.js integration for interactive visualization
- JSON data embedding with proper escaping
- Responsive design with scale switching controls
- Hover tooltips with detailed information

#### utils.c
**Purpose**: Utility functions and file operations
**Key Functions**:
- `mkdir_p()`: Recursive directory creation
- `file_exists()`: File existence checking
- `get_basename()`: Cross-platform basename extraction
- `hash_int32()`, `equal_int32()`: Hash table utilities for 32-bit keys

#### queue.c
**Purpose**: Thread-safe queue implementation for BFS operations
**Key Functions**:
- `init_queue()`: Initialize queue with dynamic sizing
- `enqueue()`, `dequeue()`: Thread-safe queue operations
- `is_empty()`, `clear_queue()`: Queue state management

**Technical Details**:
- Circular buffer implementation with automatic resizing
- Used for connected component analysis in UMI deduplication
- Lock-free operations for high-performance graph traversal

#### globals.c
**Purpose**: Global variable definitions and lookup table initialization
**Key Variables**:
- `seq2code[]`, `code2seq[][]`: DNA sequence encoding/decoding tables
- `diff2Hamming[]`: Hamming distance lookup for XOR differences
- `whitelist`, `whitelist_hash`: Barcode whitelist storage
- `feature_code_hash`: Global feature code lookup table

---

## Script Files (scripts/)

### Testing and Validation Scripts

#### test_assign.sh
**Purpose**: Integration test for assignBarcodes with realistic parameters
**Usage**: `./test_assign.sh`
**Features**:
- Tests full pipeline with NXT barcode translation
- Uses production-like parameters (limit_search=2, max_barcode_mismatches=1)
- Validates against known reference feature set

#### test_demux.sh  
**Purpose**: Comprehensive test suite for demux_fastq
**Features**:
- Tests both hash and direct search modes
- Input/output read count validation
- Per-file-set consistency checking
- Performance timing and comparison

#### run_demux.sh
**Purpose**: Production demux_fastq runner with error checking
**Features**:
- Environment variable configuration
- Input validation and error handling
- Performance benchmarking with timing
- Output verification and consistency checks

### Utility Scripts

#### downsample_fastq_directory.sh
**Purpose**: FASTQ file downsampling for testing and development
**Usage**: `./downsample_fastq_directory.sh <nReads> [dir]`
**Features**:
- Preserves gzip compression in outputs
- Handles both .fastq and .fastq.gz files
- Creates downsampled/ subdirectory with same file structure
- SIGPIPE-safe gzip handling

**Technical Details**:
- Uses `head -n` for precise line counting (4 lines per read)
- Maintains file extensions and compression
- Robust error handling with `set -euo pipefail`

#### offset_sweep.py
**Purpose**: Optimize probe detection offset parameters
**Usage**: `python3 offset_sweep.py --test_dir <dir> --offset <center> --window <range>`
**Features**:
- Samples configurable number of reads for analysis
- Tests offset ranges around expected position
- Generates tabular output with match percentages
- Supports all read types (R1, R2, R3)

**Algorithm**:
- Loads probe k-mers from barcode table
- Extracts sequences from FASTQ files
- Tests exact matches at different offsets
- Reports match rates for optimization

#### generate_colormaps.py
**Purpose**: Generate C header files for heatmap colormaps
**Features**:
- Creates plasma colormap arrays at different resolutions
- Outputs C-compatible header files
- Supports 16, 64, 256, and 1024 color levels
- Ensures perceptually uniform color scaling

### Alignment and Processing Scripts

#### runMultiAlign.sh
**Purpose**: STAR Solo alignment wrapper with multiple tool support
**Features**:
- Configurable barcode/UMI lengths via environment variables
- Support for STAR Solo, Alevin-fry, and Kallisto
- Automatic output directory organization
- Comprehensive STAR parameters for single-cell RNA-seq

**STAR Configuration**:
- CB_UMI_Simple chemistry with configurable lengths
- EmptyDrops cell filtering
- Multi-gene UMI filtering and 1MM_CR deduplication
- Complete SAM attribute output (CB, UB, CR, CY, UR, UY, GX, GN)

#### run_demux_bam.sh
**Purpose**: BAM processing example with production parameters
**Features**:
- Configurable probe offset and thread counts
- MAPQ filtering and verbose output
- Example of typical demux_bam usage patterns

### Analysis and QC Scripts

#### calculate_signal.sh
**Purpose**: Signal-to-noise analysis for feature assignment results
**Features**:
- Processes Matrix Market output files
- Calculates feature-specific signal metrics
- Generates summary statistics for QC assessment

#### count_unique_mtx_values.sh
**Purpose**: Matrix Market file analysis and validation
**Features**:
- Counts unique values in sparse matrices
- Validates Matrix Market format compliance
- Generates distribution statistics for count data

#### assign_MSK_gene.sh
**Purpose**: Memorial Sloan Kettering gene assignment workflow
**Features**:
- Institution-specific parameter sets
- Optimized for MSK sequencing protocols
- Production-ready parameter validation

---

## Key Design Patterns and Architectures

### Memory Management
- **Pool-based allocation**: Eliminates malloc/free overhead for frequent allocations
- **Block-based expansion**: Reduces memory fragmentation
- **Thread-local pools**: Minimizes contention in multi-threaded contexts
- **Aligned structures**: Optimizes cache performance

### Parallelization Strategy
- **Process-level**: Multiple samples processed simultaneously
- **Thread-level**: Producer-consumer model within samples
- **SIMD optimization**: Vectorized sequence operations
- **Lock-free data structures**: Ring buffers and atomic operations

### Error Handling and Robustness
- **Bayesian error correction**: Posterior probability-based barcode rescue
- **Connected component analysis**: UMI error correction through graph algorithms
- **Configurable stringency**: Multiple deduplication strategies
- **Comprehensive validation**: Input validation and format checking

### Performance Optimizations
- **2-bit sequence encoding**: 4x memory reduction for DNA sequences
- **Direct vs. hash lookup**: Automatic algorithm selection based on data size
- **Memory-mapped I/O**: Efficient file handling for large datasets
- **Lazy initialization**: On-demand resource allocation

This architecture provides a scalable, high-performance solution for single-cell feature assignment with extensive configurability and robust error handling.
