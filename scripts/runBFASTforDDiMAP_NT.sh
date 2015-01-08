#!/bin/bash
#this program expects Fastq and .fa file to exist

cd $myWorkingDir
pwd

bfast fasta2brg -f $myFasta -A 0

bfast index -f $myFasta -m 1111111111111111111111 -w 14 -i 1 -A 0 -n 8
bfast index -f $myFasta -m 111110100111110011111111111 -w 14 -i 2 -A 0 -n 8
bfast index -f $myFasta -m 10111111011001100011111000111111 -w 14 -i 3 -A 0 -n 8
bfast index -f $myFasta -m 1111111100101111000001100011111011 -w 14 -i 4 -A 0 -n 8
bfast index -f $myFasta -m 111111110001111110011111111 -w 14 -i 5 -A 0 -n 8
bfast index -f $myFasta -m 11111011010011000011000110011111111 -w 14 -i 6 -A 0 -n 8
bfast index -f $myFasta -m 1111111111110011101111111 -w 14 -i 7 -A 0 -n 8
bfast index -f $myFasta -m 111011000011111111001111011111 -w 14 -i 8 -A 0 -n 8
bfast index -f $myFasta -m 1110110001011010011100101111101111 -w 14 -i 9 -A 0 -n 8
bfast index -f $myFasta -m 111111001000110001011100110001100011111 -w 14 -i 10 -A 0 -n 8

export myBmf=$myPatient.bmf

bfast match -f $myFasta -A 0 -i 1-10 -k 18 -K 100000 -w 0 -t -n 8 -Q 100000 -l -r $myFastq > $myBmf

export myBaf=$myPatient.baf

bfast localalign -f $myFasta -m $myBmf -A 0 -n 8 -U -q 20 -Q 100000 -t > $myBaf

#  bmf file no longer needed
rm $myBmf
bfast postprocess -f $myFasta -i $myBaf -o $myAligned -O 1 -a 3 -z -n 8 -q 20 -Q 100000 -t > $myAligned.sam

#  baf, bif, brg files no longer needed
rm $myBaf
rm $myAligned*.bif
rm $myAligned*.brg
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
