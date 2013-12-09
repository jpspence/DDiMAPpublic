High Performance DDiMAP
=======================

DDiMAP in C and CUDA. 

### Running the code
Compiled code is placed in bin/.
```bash
make 
./bin/DDiMAP
```

### Installing 
```bash
make install
make
```



### Available make commands
```bash

make            # Compiles all the source code, generate executable
make install    # This will download all the dependencies
make test       # Runs all the unit tests for the project
make clean      # Removes everything

```

### Submodules 

https://github.com/pezmaster31/bamtools     # Tools for accessing BAM files in C

https://github.com/danfis/cu                # Unit testing in C | http://cu.danfis.cz/

https://github.com/boostorg/program_options # Library used to parse command line arguments.

### Resources for running on EC2
http://vasir.net/blog/opencl/installing-cuda-opencl-pyopencl-on-aws-ec2

