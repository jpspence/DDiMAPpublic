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

	// ------------------------------------------------------------------------
	// Parameter Parsing
	// ------------------------------------------------------------------------
	static struct option long_options[] = {
			{"file", 	0, 0, 'f'},
			{NULL, 		0, NULL, 0}
	};

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:", long_options, &option_index)) != -1) {

		switch (c) {
		case 'f':
			printf ("Parsing file :  %s \n",optarg);
			file = optarg;
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
	printf ("It took me %lu ticks (%f seconds) to read %d | %d reads from BAM file.\n",
			t, ((float)t)/CLOCKS_PER_SEC, unique, total);



	t = clock();
	sequential();
	int verified = printFasta();
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to verify %d | %d.\n",
			t, ((float)t)/CLOCKS_PER_SEC , verified, unique);


	t = clock();
	verified = iterate(buildHistograms);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to build the histogram.\n",
			t, ((float)t)/CLOCKS_PER_SEC);

	printHistograms();

	callSNVs();

	gnuplot_ctrl    *   h1;
	h1 = gnuplot_init() ;
	gnuplot_resetplot(h1) ;
	gnuplot_cmd(h1, "set xrange [0:1200]");
	gnuplot_cmd(h1, "set yrange [0.0001:0.9999]");
	gnuplot_cmd(h1, "set multiplot title \"GC Content for different genes\" layout 3,4");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Bcl2_NCBICG.txt\"     with points ls 1 title \"Bcl2\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Bcl6_NCBICG.txt\"     with points ls 1 title \"Bcl6\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/CD83_NCBICG.txt\"     with points ls 1 title \"CD83\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Emu_NCBICG.txt\"      with points ls 1 title \"Emu\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/J6-J4_NCBICG.txt\"    with points ls 1 title \"j6-j4\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/KS_NCBICG.txt\"       with points ls 1 title \"KSCG\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Pax5_NCBICG.txt\"     with points ls 1 title \"Pax5\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Pim1_NCBICG.txt\"     with points ls 1 title \"Pim1\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/RhoH_NCBICG.txt\"     with points ls 1 title \"RhoH\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Vh1-18_NCBICG.txt\" with points ls 1 title \"vh1\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/myc-1_NCBICG.txt\"    with points ls 1 title \"myc-1\"");
	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/myc-2_NCBICG.txt\" with points ls 1 title \"myc-2\"");
	gnuplot_cmd(h1, "unset multiplot") ;

	cin.get();

	// ------------------------------------------------------------------------
	// End.
	// ------------------------------------------------------------------------

	return EXIT_SUCCESS;

}
