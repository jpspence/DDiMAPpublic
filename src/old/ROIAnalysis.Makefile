C = gcc
CCFLAGS = -g -O3
LDFLAGS = -g -O3
OUTPUTDIR = ./bin/

all: ROIAnalysis ../samtools/libbam.a ../samtools/samtools ../bcftools/bcftools \
		 ex1.glf ex1.pileup.gz ex1.bam.bai ex1f-rmduppe.bam ex1f-rmdupse.bam ex1.glfview.gz ex1.bcf calDepth
     @echo; echo \# You can now launch the viewer with: \'samtools tview ex1.bam ex1.fa\'; echo;

ROIAnalysis: ../samtools/libbam.a ROIAnalysis.o
	$(CC) $(LDFLAGS) -o $(OUTPUTDIR)ROIAnalysis $(OUTPUTDIR)ROIAnalysis.o -L.. -../samtools/lbam -lm -lz

ROIAnalysis.o: ROIAnalysis.c bin
	$(CC) $(CCFLAGS) -c ROIAnalysis.c -o $(OUTPUTDIR)ROIAnalysis.o

bin:
	mkdir -p $(OUTPUTDIR)

clean:
	-rm  -rf $(OUTPUTDIR)

../samtools/samtools:
	(cd ../samtools; make samtools)

../samtools/libbam.a:
	(cd ../samtools; make libbam.a)

clean:
	rm -fr *.bam *.bai *.glf* *.fai *.pileup* *~ calDepth *.dSYM ex1*.rg ex1.bcf
