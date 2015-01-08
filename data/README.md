High Performance DDiMAP Test Examples
=====================================

Small data sets for testing code installation and iterative script testing

### Running a test

# Open a terminal or other command window

# Make sure DDiMAP, samtools, bamtools, and the mapper(s) you are testing are on your path (your locations will vary)

which DDiMAP samtools bamtools bfast bwa cushaw3 novoalign gmapper  #  returns a path to each executable found in your PATH

export PATH=$PATH:$HOME/DDiMAP/bin          #  required to use the DDiMAP executable and the DDiMapIterate.py script
export PATH=$PATH:$HOME/samtools-0.1.19     #  required to use the DDiMapIterate.py script

export PATH=$PATH:$HOME/bwa-0.7.10          #  required only if bwa mem is to be used in mapping and not already in path
export PATH=$PATH:$HOME/SHRiMP_2_2_3/bin    #  required only if shrimp2 (gmapper) is to be used in mapping
export PATH=$PATH:$HOME/cushaw3-v3.0.2      #  required only if cushaw3 is to be used in mapping
export PATH=$PATH:$HOME/novocraft           #  required only if novoalign is to be used in mapping

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/bamtools  #  pointing it to the bamtools library if needed, this is the make install location

# cd into /path/to/DDiMAP/data (your installation location will vary)

cd $HOME/DDiMAP/data

# cd into a test directory (CStest or IterativeTest or PairedTest) and run test using python to run script

python ../../DDiMapIterate.py

# different test cases are in different cfg files, use -c to specify file to use (./DDiMAP.cfg is default)

python ../../DDiMapIterate.py -c ./DDiMAP.cfg

