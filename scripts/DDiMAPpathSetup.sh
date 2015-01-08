#!/bin/bash
#
#  edit to reflect the paths to the executables that you want to use with DDiMapIterate.py
#
export PATH=$PATH:$HOME/DDiMAP/bin          #  required to use the DDiMAP executable and the DDiMapIterate.py script
export PATH=$PATH:$HOME/samtools-0.1.19     #  required to use the DDiMapIterate.py script
export PATH=$PATH:$HOME/bwa-0.7.10          #  required only if bwa mem is to be used in mapping, may comment out using initial # if not used
export PATH=$PATH:$HOME/SHRiMP_2_2_3/bin    #  required only if shrimp2 is to be used in mapping
export PATH=$PATH:$HOME/cushaw3-v3.0.2      #  required only if cushaw3 is to be used in mapping
export PATH=$PATH:$HOME/novocraft           #  required only if novoalign is to be used in mapping

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/bamtools  #  pointing it to the bamtools library is needed, this is the make install location



