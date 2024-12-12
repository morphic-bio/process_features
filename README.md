## Introduction

Biodepot assignBarcodes is a fast parallel utility (**1 Million reads per/second on a laptop**) for targeted sequencing) to quickly assign a set of feature barcodes that have been sequenced to a set sequence barcodes for single cell analysis. It is a much faster more powerful and open-source alternative to the custom feature option in CellRanger for this purpose. Unlike CellRanger, assignBarcodes performs an exhaustive search which allows it to find feature barcodes in ATAC seq and expressed in RNA-seq fastq files, and increases the coverage in targeted sequencing. The same sequence barcode error correction methodologies employed by CellRanger are used by default and can be customized. Error correction and fuzzy matching options are available for feature sequence matching. When available, counts are deduped using UMIs and the user can specify different options to handle how counts are collapsed depending on whether targeted sequencing is being used. UMI sequencing errors are anticipated and handled by creating connected cliques of similar UMIs before deduping. 

## Quick start
### Installation

On Ubuntu:
Clone and cd into this repository
```
sudo apt-get update
sudo apt-get -y install build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libgomp1 libglib2.0-dev libcairo2-dev
cd process_features
make
sudo cp assignBarcodes /usr/local/bin
```
If the heatmap png is not needed, you can leave out libcairo2-dev and compile with
```
make NO_HEATMAP=1
```

#### Docker
Alternatively you can use the Docker container biodepot/process_features:latest or build it using the Dockerfile provided

#### Usage
assignBarcodes -w whitelist.txt -f features.txt -b 16 -m 2 -s 10 -d output/ --barcode_fastqs barcodes.fastq --forward_fastqs forward.fastq --reverse_fastqs reverse.fastq


## Flags and Descriptions

### Input and Output
- `-w, --whitelist [filename]`: Specify the whitelist file for barcodes.
- `-f, --featurelist [filename]`: Provide the feature list file.
- `-d, --directory [path]`: Set the working directory for output files.

### Barcode and Feature Settings
- `-b, --barcode_length [int]`: Specify the length of the barcodes.
- `-u, --umi_length [int]`: Set the length of the UMI.
- `-o, --feature_constant_offset [int]`: Constant offset for feature processing.
- `-B, --barcode_constant_offset [int]`: Constant offset for barcode processing.
- `--barcode_fastqs [string]`: Provide barcode FASTQ files.
- `--forward_fastqs [string]`: Provide forward FASTQ files.
- `--reverse_fastqs [string]`: Provide reverse FASTQ files.
- `--barcode_fastq_pattern [string]`: Set the pattern for barcode FASTQ files.
- `--forward_fastq_pattern [string]`: Set the pattern for forward FASTQ files.
- `--reverse_fastq_pattern [string]`: Set the pattern for reverse FASTQ files.

### Error Handling and Thresholds
- `-m, --maxHammingDistance [int]`: Maximum Hamming distance allowed for barcodes.
- `-s, --stringency [int]`: Stringency level for barcode matching.
- `-i, --min_counts [int]`: Minimum counts threshold.
- `-M, --min_posterior [float]`: Minimum posterior threshold.
- `--max_barcode_mismatches [int]`: Maximum mismatches allowed for barcodes.

### Performance and Threads
- `-t, --threads [int]`: Maximum number of concurrent processes (maximum sets of files that can be processed at once)
- `-T, --available_threads [int]`: Maximum available threads for processing.
- `-S, --search_threads [int]`: Threads allocated for searching barcodes.
- `-c, --consumer_threads_per_set [int]`: Consumer threads per data set.
- `-R, --read_buffer_lines [int]`: Number of lines to buffer when reading files.
- `-L, --average_read_length [int]`: Average length of reads.

### Boolean Flags
- `-v, --debug`: Enable debug mode.
- `-a, --as_named`: Process files in named mode.
- `-r, --reverse_complement_whitelist`: Reverse complement the whitelist.
- `-k, --keep_existing`: Preserve existing outputs.
- `-P, --parallel_by_file`: Enable parallel processing by file.

### Other Settings
- `--feature_n [int]`: Maximum feature count.
- `--barcode_n [int]`: Maximum barcode count.

## Example Usage

```bash
./assignBarcodes -v -T 32 -R 1024 --barcode_fastq_pattern _R1_ --forward_fastq_pattern _R2_ --reverse_fastq_pattern _R3_ -w  /mnt/pikachu/processing/genome/10x_version3_whitelist.txt -o 26 -m 5 -t 7 -s 1 -i 0 -u 12  --max_barcode_mismatches 3  -f /mnt/pikachu/project_14361_feature_ref_v2_202405_alt.csv -d msk_test_heatmap  -d msk_targeted_2 /mnt/pikachu/morphic-bio-processing/morphic-msk/pooled_scRNA_seq_202401/larry-barcodes/A/
```


## Search methodology
#### Initial fixed position search
For targeted sequencing, most of the of read sequences will have a constant start sequence of fixed length. assignBarcodes attempts to make a match here first. If no match is found, then it does a more expensive exhaustive search.
#### Exhaustive search
The exhaustive search checks the entire read against all the feature barcodes at all possible starting positions in the read. For ATAC-seq the search is done in both orientations. We use a novel method that converts the query and match sequences to bitcodes. We uses bitwise ops and a lookup table for hamming distance evaluation of 4 basepairs chunks with a bitops and lookup and can be vectorized by the compiler for even greater speedup. Additionally, the search is broken down into four independent subsearches which are performed in parallel for a 16x speedup over the simple Hamming search.
## Error correction
### Sequence barcodes
The error correction handles Ns (unknown base pairs) and sequencing errors. To take into account sequencing errors, a barcode can be at most 1 base pair different from a single valid barcode and then it will be assigned to that barcode. If there are multiple barcodes, then we look at the quality scores and the number of barcodes variants observed and find the most likely match for the barcode based on the posterior probability. This is described in the Cell Ranger documentation. 
To handle N's the user specifies a maximum number of Ns (at compile time - default is 1) that are tolerated. All the possible base pairs are substituted for an N and then compared to see if a unique barcode is found. 
### Feature barcodes
To handle sequencing errors, the user specifies a maximum Hamming distance. If a sequence matches a feature barcode within the Hamming distance and uniquely to a sequence with a minimum distance then it is assigned to that feature barcode. For N's up to a maximum specified by the user, all possible variations are generated for the N's and checked against the possible sequences. If there is a unique best match (minimum Hamming distance) that is less or equal to the maximum Hamming distance then it is assigned to that feature barcode. Assignments are tentative, pending the completion of the comprehensive search (unless there is an exact match). If there is no exact match, the comprehensive search attempts to find a better match.

## UMI de-duplication

#### Introduction
Unique Molecular Identifiers are random sequences that are included with barcodes that are used to identify groups of reads that are duplications of the same sequence. De-duping counts using gives a more accurate reflection of the relative abundance of the originating sequences. 
The de-duping algorithm depends on the methodology used.
In scRNA-seq and bulk RNA-seq, counts assigned to a sequence barcode-umi are collapsed only if they map to the same position. This is a relatively rare event that occurs only due to the large number of reads. For CRISPR-targeted sequencing, there is a much lower number of possible mappings and it is very common for a barcode-umi to map to multiple feature sequences and requires a more complicated strategy. assignBarcodes gathers counts all the sequences associated with a barcode-umi, and chooses a user option to handle how the counts should be de-duped. 
#### Aggregation into connected components
Sequencing errors of UMIs can occur and single sequence errors in UMIs are much more likely than having two random UMIs that differ by 1 base pair. To account for this, assignBarcodes aggregates barcode-umi sets that differ by 1 base pair in the UMI (connected-component). 

Once the connected component is formed. The counts are aggregated based on two variables, stringency and minimum_counts that are provided by the user using the -s and -i flags.

#### De-duping strategy options
#### stringency=0
When this is the case, we use the RNA-seq strategy. Any feature counts are reduced to 1 if there are at least min_counts. This is added to the feature counts associated with the sequence barcode. This is the only case where more than one feature is assigned a count from a barcode-umi set.
#### stringency=1-999
When this is the case, we find the feature with the highest count. If there is a tie then we assign no feature counts to the barcode. Otherwise, we look at the feature with the maximum number of counts assigned to it. If this number of counts is greater than total_counts *( counts/1000) and total_counts > min_counts then we assign a feature count to the barcode.
#### stringency >=1000
We assign a feature count only if there is one feature with counts and the counts are greater than min_counts.

### Parallelization

There are 3 levels of parallelization. The user can completely control the number of threads used in each parallelization process or can just assign the number of threads and allow assignBarcodes to apportion them out optimally. The different levels of parallelization used are:
1.  A producer-consumer model that allows the reading of fastq files to take place simultaneously with the processing
2. Parallelization of the hamming search
3. Forking off separate processes for multiple samples

#### Fastq input and processing
This is done using pthreads. The default strategy is to have one thread read the barcode file, the forward file and if it exists the reverse read file. The user also has the option to use a thread for each fastq file. Though this has more overhead, it can be faster when there is are only a few samples that are split into many files and/or there are many processors available with a fast SSD. One or more consumer threads pull reads from the buffer and processes them while the fastq reader(s) populate them from the files.

#### Hamming search parallelization

The search function proceeds using 4 independent search frames. OpenMP is used to assign up to 4 threads to parallelize and increase the search speed.
#### Forking off separate processes
To handle multiple samples, assignBarcodes can fork off a separate process. This is considerably faster than using OpenMP with little extra overhead as the memory footprint is very small (often less than 1GB per process - closer to 100 MB for ATAC and scRNA-seq files as most of their reads do not have feature barcodes).
### Feature and sequence file formats

Sequence barcodes are provided one barcode per line. Standard sequence barcode whitelists from 10x work fine.

Feature barcodes should be provided as comma separated files with a header line. The header line must contain a 'sequence' field and a 'name' field. The other fields are ignored.

### Output format
The assignments are outputted in matrix market format which essentially has 3 files listing the barcodes, features and count matrix.

### QC files
#### Run stats
In the output directory some statistics from the run are kept in "stats.txt"
```
Total feature counts 31993392
Total deduped feature counts 6765455
Total unique barcode UMIs 9148050
Total whitelisted barcodes 259694
Total_unmatched_reads 7236844
Percentage reads assigned to barcode 81.5529
```
#### Matched sequences
Each assigned feature sequence is listed in the "feature_sequences.txt" file.
```
Feature Index Sequence Hamming Distance Counts Feature Name  
  5 CCAGTGCAGTCACTCGACACCTGAGATCGTGACCAGACCA  0 1267629 17_ARX
  5 nCAGTGCAGTCACTCGACACCTGAGATCGTGACCAGACCA  0       9 17_ARX
  5 CnAGTGCAGTCACTCGACACCTGAGATCGTGACCAGACCA  0       8 17_ARX
  5 nnAGTGCAGTCACTCGACACCTGAGATCGTGACCAGACCA  0       6 17_ARX
  5 nCnGTGCAGTCACTCGACACCTGAGATCGTGACCAGACCA  0       5 17_ARX
  5 CCAGTGCAGTCnCTCGACACCTGAGATCGTGACCAGACCA  0       1 17_ARX
  5 CCAGTGCAGTCnCTCGtCACCTGAGATCGTGtCCAGACCA  0       1 17_ARX
  5 gCcagtgtagtcactcgacaCctgagatcgtgaCcagaCc  1   10115 17_ARX
  5 gCcagtgcagtcactcgacaCctgagatcgtgaCcagaCc  1    3911 17_ARX
...
```
Each of the matched sequences is displayed under the feature index that they are matched to. Currently, the optimal alignment is not displayed - they are just shown from the start and only up to the length of the feature sequence - but this may change in the future. 
#### Feature frequency histogram heatmap
After deduplication, we can make a histogram of the number of features assigned to cells. The perfect scenario would be nearly all cells have 1 feature assigned to them with a small distribution representing multiplets and another small distribution from empty cells (ambient counts that are assigned a cell barcode). Underneath is a heatmap where each row is the same histogram, but only showing cells that have that feature in it. The values of the histograms are space separated in feature_histograms.txt. There is also another file with the coexpression matrix, i.e. how many times feature i is seen with feature j. This is in feature_coexpression.txt. Currently this is not plotted.

###### Heatmap example
![Heatmap example](./graphics/heatmap.png)


