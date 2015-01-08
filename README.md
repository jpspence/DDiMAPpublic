High Performance DDiMAP
=======================
DDiMAP in C

### Prerequisites
You'll need the following tools : `git` `gcc` `make` `cmake`

### Installing the code
```bash
#This currently doesn't run on Windows, I have no plans on supporting it either.
make setup
make
```

### Running the code

```bash

# Run DDiMAP on a fasta and bam file
./bin/DDiMAP -f /location/of/fasta.fa  -b /location/of/bam.bam

# View usage details
./bin/DDiMAP --help
usage : DDiMAP [-f <fasta> -b <bam> <args>] [--help] [--test]

Basic Parameters:
   --bam              | -b   This specifies the path to the bam file
   --fasta            | -f   This specifies the path to the fasta file
   --keepID           | -k   Keep reads that have both an insert and delete in CIGAR string
   --roa-size         | -r   Number of base pairs for a Region of Analysis    (default : 34)
   --min-count-thresh | -c   Minimum count of reads to see in each direction (default : 2)
   --min-freq-thresh  | -p   Minimum frequency of reads to see in each direction (default : 0.00075)

Frag Making Parameters
   --frag-threshold   | -a   Minimum verified word frequency threshold for use in fragment assembly (default : 0.01)
   --nv-threshold     | -n   Minimum non-verified word coverage for use in fragment assembly (default : 0.1)

SNV Call Reason Parameters
   --snv-verified     | -s   Minimum observed variant frequency in verified words for a CallReason 2 SNV (default : 0.003)
   --snv-total        | -l   Minimum observed variant frequency in all words for a CallReason 3 SNV (default : 0.1)

Output Parameters
   --output           | -o   Directory to store output (default : ./output/)
   --dictionary-level | -d   Dictionary verbosity : 0 = fwd/rev counts only | 1 = add full read in/del counts | 2 = add frag counts (default : 1)

### Available make commands
```bash

make            # Compiles all the source code, generate executable
make setup      # This will download all the dependencies
make clean      # Removes DDiMAP executable and *.o files

```

### Submodules used by this repository

https://github.com/pezmaster31/bamtools     # Tools for accessing BAM files in C
https://github.com/lh3/readfq               # Tools for reading fasta and fastq files
