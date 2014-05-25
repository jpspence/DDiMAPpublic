High Performance DDiMAP
=======================
DDiMAP in C and CUDA.

### Prerequisites
You'll need the following tools : `git` `gcc` `make` `cmake` `cuda`

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

# Run DDiMAP's test suite
./bin/DDiMAP --test

# View usage details
./bin/DDiMAP --help
usage : DDiMAP [-f <fasta> -b <bam> <args>] [--help] [--test]

Basic Parameters:
   --bam              | -b   This specifies the path to the bam file
   --fasta            | -f   This specifies the path to the fasta file
   --keepID           | -k   Keep reads that have both an insert and delete in CIGAR string
   --verify-threshold | -v   Minimum number of reads to see in each direction (default : 2)
   --roa-size         | -r   Number of base pairs for a Region of Analysis    (default : 34)

Frag Making Parameters
   --ppm              | -p   Minimum level of reads to consider for DDiMAP    (default : 0.00075)
   --frag-threshold   | -a   Minimum verified coverage required to be considered for frags (default : 0.01)
   --nv-threshold     | -n   Minimum non-verified coverage required to be considered for frags (default : 0.1)

SNV Calling Parameters
   --snv-verified     | -s   Minimum level of nucleotide variation in verified words to call an SNV (default : 0.003)
   --snv-total        | -l   Minimum level of nucleotide variation in total to call an SNV (default : 0.1)

Output Parameters
   --output           | -o   Directory to store output (default : ./output/)
   --dictionary-level | -d   Dictionary verbosity : 0 = fwd/rev counts | 1 = in/del data | 2 = frag mappings (default : 0)

Testing Parameters
   --test             | -t   Run the test suite

Future Parameters (works in progress):
   --length-of-snv-ref| -l   Number of base pairs you'd like to see in SNV
```


### Available make commands
```bash

make            # Compiles all the source code, generate executable
make setup      # This will download all the dependencies
make clean      # Removes everything

```

### Submodules used by this repository

https://github.com/pezmaster31/bamtools     # Tools for accessing BAM files in C

https://github.com/danfis/cu                # Unit testing in C | http://cu.danfis.cz/


### Resources for running on EC2
http://vasir.net/blog/opencl/installing-cuda-opencl-pyopencl-on-aws-ec2

