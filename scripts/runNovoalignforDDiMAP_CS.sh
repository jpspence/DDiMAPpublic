#!/bin/bash
#this program expects Fastq and .fa file to exist

cd $myWorkingDir
pwd

/home/jspence/novocraft/novoindex -c ${myAligned}.nix $myFasta
/home/jspence/novocraft/novoalignCS -r R -o SAM -d ${myAligned}.nix -f $myFastq > ${myAligned}.sam
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
samtools view -b -S -o $bam_file $sam_file
#  sam file no longer needed
rm $sam_file
#  sort sorts them along chromsomes to make index work
#echo "Making sorted bam file"
#samtools sort $bam_file $bamsorted_prefix
#  bam file no longer needed
#rm $bam_file
#  now make the index
#echo "Indexing sorted bam file"
#samtools index $bamsorted_file
