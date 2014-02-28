//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP in C
//============================================================================

#include <getopt.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "DDiMAP-lib.h"
#include "../include/gnuplot_i/gnuplot_i.h"

#define SLEEP_LGTH  1

// Default file.

//string file  = "data/128test_Gen1_example_sorted.bam";
//char  *fasta = "data/128test_Gen1_example.fa";

string file  = "/Dropbox/Google Drive/DataExchangeUR/128_Gen7_BSBSB_VhJ_FragHigh_sorted.bam";
char  *fasta = "/Dropbox/Google Drive/DataExchangeUR/128_Gen7_BSBSB_VhJ_FragHigh.fa";

int main (int argc, char **argv)
{
	int c;
	clock_t t;
	int total = 0, unique = 0;

	// Default values
	// FLAGS
	// TODO: -dropID
	// TODO: -uniqueThreshold
	int VERIFY_THRESHOLD  =  2;
	double PPM = 0.00075;
	double FRAG_THRESHOLD = .01;
	double SNV_VERIFIED_THRESHOLD = .003;
	double SNV_TOTAL_THRESHOLD = .1;

	// ------------------------------------------------------------------------
	// Parameter Parsing
	// ------------------------------------------------------------------------
	static struct option long_options[] = {
			{"file", 	0, 0, 'f'},
			{"verify-threshold", 	0, 0, 'v'},
			{"ppm", 	0, 0, 'p'},
			{"fragment-threshold", 	0, 0, 'a'},
			{"frequency-threshold", 	0, 0, 'b'},
			{"snv-threshold", 	0, 0, 'c'},
			{NULL, 		0, NULL, 0}
	};

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:a:b:c:v:p:", long_options, &option_index)) != -1) {

		switch (c) {
		case 'f':
			printf ("Parsing file :  %s \n",optarg);
			file = optarg;
			break;
		case 'a':
			printf ("Applying threshold of :  %s \n",optarg);
			break;
		case 'b':
			printf ("Applying threshold of :  %s \n",optarg);
			break;
		case 'c':
			printf ("Applying threshold of :  %s \n",optarg);
			break;
		case 'v':
			VERIFY_THRESHOLD = atoi(optarg);
			printf ("Applying threshold of :  %d \n",VERIFY_THRESHOLD);
			break;
		case 'p':
			PPM = atof(optarg);
			printf ("Applying threshold of :  %f \n",PPM);
			break;
		default:
			printf ("?? getopt returned character code 0%o ??\n", c);
		}
	}
	if (optind < argc) {
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		printf ("\n");
	}


	// ------------------------------------------------------------------------
	// DDiMAP
	// ------------------------------------------------------------------------

	t = clock();
	unique = readFile(file, fasta, 34, buildRead);
	t = clock() - t;
//	printf ("It took me %lu ticks (%f seconds) to read %d | %d reads from BAM file.\n",
//			t, ((float)t)/CLOCKS_PER_SEC, unique, total);



	t = clock();
	sequential(VERIFY_THRESHOLD, PPM, FRAG_THRESHOLD);
	int verified = printFasta();
	t = clock() - t;
//	printf ("It took me %lu ticks (%f seconds) to verify %d | %d.\n",
//			t, ((float)t)/CLOCKS_PER_SEC , verified, unique);


	t = clock();
	verified = iterate(buildHistograms);
	t = clock() - t;
//	printf ("It took me %lu ticks (%f seconds) to build the histogram.\n",
//			t, ((float)t)/CLOCKS_PER_SEC);

	printHistograms();

	iterate(count);

	callSNVs(SNV_VERIFIED_THRESHOLD, SNV_TOTAL_THRESHOLD);

	gnuplot_ctrl    *   h1;
	h1 = gnuplot_init() ;
	gnuplot_resetplot(h1) ;
	gnuplot_cmd(h1, "set xrange [0:1200]");
	gnuplot_cmd(h1, "set yrange [0.0001:0.9999]");
	gnuplot_cmd(h1, "set multiplot title \"GC Content for different genes\" layout 3,4");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Bcl2CG.txt\"     with points ls 1 title \"Bcl2\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Bcl6CG.txt\"     with points ls 1 title \"Bcl6\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/CD83CG.txt\"     with points ls 1 title \"CD83\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/EmuCG.txt\"      with points ls 1 title \"Emu\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/J6-J4CG.txt\"    with points ls 1 title \"j6-j4\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/KSCG.txt\"       with points ls 1 title \"KSCG\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Pax5CG.txt\"     with points ls 1 title \"Pax5\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Pim1CG.txt\"     with points ls 1 title \"Pim1\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/RhoHCG.txt\"     with points ls 1 title \"RhoH\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Vh1-18CG.txt\" with points ls 1 title \"vh1\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/myc-1CG.txt\"    with points ls 1 title \"myc-1\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/myc-2CG.txt\" with points ls 1 title \"myc-2\"");
	gnuplot_cmd(h1, "unset multiplot") ;

	cin.get();

	// ------------------------------------------------------------------------
	// End.
	// ------------------------------------------------------------------------

	return EXIT_SUCCESS;

}
