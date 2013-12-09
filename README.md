High Performance DDiMAP
=======================
DDiMAP in C and CUDA. 

### Prerequisites
As much of this code has been modularized to run on vanilla installs, assuming
that you have the following `git` `gcc` `make` `cmake` `cuda`

I have tried my best to have a self contained repo that should be able to be 
brought up as quickly as possible.

### Installing the code
This currently doesn't run on Windows, and I have no plans on support it either.
```bash
make setup 
make
```

### Running the code
All of the compiled code is placed in `bin`
```bash
make 
# Optional -f Flag allows you to provide a file to process
./bin/DDiMAP    # This runs the serial code.
./bin/DDiMAPGPU # This runs the GPU


```


### Available make commands
```bash

make            # Compiles all the source code, generate executable
make setup      # This will download all the dependencies
make test       # TODO: Runs all the unit tests for the project
make clean      # Removes everything

```

### Submodules used by this repository

https://github.com/pezmaster31/bamtools     # Tools for accessing BAM files in C

https://github.com/danfis/cu                # Unit testing in C | http://cu.danfis.cz/


### Resources for running on EC2
http://vasir.net/blog/opencl/installing-cuda-opencl-pyopencl-on-aws-ec2

