#!/bin/bash
#this program expects Fastq and .fa file to exist

cd $myWorkingDir
pwd

export shrimpPath=/home/jspence/SHRiMP_2_2_3/bin/gmapper-cs

# runs shrimp to get sam file using 8 processors and fastq input...
#     Add a log file as a sink for the 2 output like 2>$myAligned.log if you would like (Matlab can be used to capture all stuff and log it)
#         use parameters for
#		8 processors (-N 8)
#		fastq input (-Q)
#		output only one per read (-o 1 --strata)
#		all contigs are present (--all-contigs) - this results in smaller sam file as no split genome merge is needed
#		 
$shrimpPath -N 8 -Q -o 1 --strata --all-contigs $myFastq $myFasta > ${myAligned}.sam

#
#  define files to read or create to make an indexed bam file from a sam file
#
#
#  change the lines below to match how your files are named for this reference file and sam file combo
#
export sam_file=${myAligned}.sam
#
#  change the lines below to match how you would like the output file names to look 
# 
export bam_file=${myAligned}.bam
export bamsorted_prefix=${myAligned}_sorted
export bamsorted_file=$bamsorted_prefix.bam
#
#  view converts sam to bam
#
echo "Making unsorted bam file"
/home/jan/samtools-0.1.19/samtools view -b -S -o $bam_file $sam_file
#  sam file no longer needed
rm $sam_file
#  sort sorts them along chromsomes to make index work
#echo "Making sorted bam file"
#/home/jan/samtools-0.1.19/samtools sort $bam_file $bamsorted_prefix
#  bam file no longer needed
#rm $bam_file
#  now make the index
#echo "Indexing sorted bam file"
#/home/jan/samtools-0.1.19/samtools index $bamsorted_file
