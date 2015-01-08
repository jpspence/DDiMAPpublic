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
#include <sys/types.h>
#include <sys/stat.h>

#include "DDiMAP-lib.h"
//#include "../include/gnuplot_i/gnuplot_i.h"

#define SLEEP_LGTH  1

// Default file.
string file  = "/Users/androwis/Dropbox/documents/DataExchangeUR/128_Gen7_CUSHAW2_VhJ_sorted.bam";
string fasta = "/Users/androwis/Dropbox/documents/DataExchangeUR/128_Gen7_CUSHAW2_VhJ.fa";
string output = "./output/";

// Default values
int ROA_SIZE 	= 34;
int  MIN_COUNT_THRESH  =  2;
bool DROPID = true;

// Frag making thresholds
double MIN_FREQ_THRESH = 0.00075;
double FRAG_THRESHOLD = .01;
double NON_VERIFIED_THRESHOLD = .1;

// SNV calling thresholds
double SNV_VERIFIED_THRESHOLD = .003;
double SNV_TOTAL_THRESHOLD = .1;


// Output Parameters
int DICTIONARY_LEVEL = 1;

void usage()
{
	cout << "usage : DDiMAP [-f <fasta> -b <bam> <args>] [--help]" << endl;

	cout << endl << "Basic Parameters:" << endl;
	cout << "   --bam              | -b   This specifies the path to the bam file" << endl;
	cout << "   --fasta            | -f   This specifies the path to the fasta file" << endl;
	cout << "   --keepID           | -k   Keep reads that have both an insert and delete in CIGAR string" << endl;
	cout << "   --roa-size         | -r   Number of base pairs for a Region of Analysis    (default : " << ROA_SIZE << ")" << endl;
	cout << "   --min-count-thresh | -c   Minimum count of reads to see in each direction (default : " << MIN_COUNT_THRESH<< ")" << endl;
	cout << "   --min-freq-thresh  | -p   Minimum frequency of reads to see in each direction (default : "<< MIN_FREQ_THRESH << ")" << endl;

	cout << endl << "Frag Making Parameters" << endl;
	cout << "   --frag-threshold   | -a   Minimum verified word frequency threshold for use in fragment assembly (default : "<< FRAG_THRESHOLD <<")" << endl;
	cout << "   --nv-threshold     | -n   Minimum non-verified word coverage for use in fragment assembly (default : "<< NON_VERIFIED_THRESHOLD <<")" << endl;

	cout << endl << "SNV Call Reason Parameters" << endl;
	cout << "   --snv-verified     | -s   Minimum observed variant frequency in verified words for a CallReason 2 SNV (default : " << SNV_VERIFIED_THRESHOLD<<")" << endl;
	cout << "   --snv-total        | -l   Minimum observed variant frequency in all words for a CallReason 3 SNV (default : " << SNV_TOTAL_THRESHOLD<< ")" << endl;

	cout << endl << "Output Parameters" << endl;
	cout << "   --output           | -o   Directory to store output (default : "<< output <<")" << endl;
	cout << "   --dictionary-level | -d   Dictionary verbosity : 0 = fwd/rev counts only | 1 = add full read in/del counts | 2 = add frag counts (default : "<<DICTIONARY_LEVEL<<")" << endl;

	cout <<endl;

}

//void test()
//{
//   not doing anything so comment out here and at place called from
//}

int main (int argc, char **argv)
{
	int c;
	clock_t t;
	int total = 0, unique = 0;

	// ------------------------------------------------------------------------
	// Parameter Parsing
	// ------------------------------------------------------------------------
	static struct option long_options[] = {
			{"bam", 0,0, 'b'},
			{"dictionary_level", 0,0, 'd'},
			{"fasta", 	0, 0, 'f'},
			{"keepID", 0,0, 'k'},
			{"min-count-thresh", 	0, 0, 'c'},
			{"roa-size", 	0, 0, 'r'},
			{"output", 0,0, 'o'},

			{"min-freq-thresh", 	0, 0, 'p'},
			{"frag-threshold", 	0, 0, 'a'},
			{"nv-threshold", 	0, 0, 'n'},

			{"snv-verified", 	0, 0, 's'},
			{"snv-total", 	0, 0, 'l'},

			{"help", 	0, 0, 'h'},
			{NULL, 		0, NULL, 0}
	};

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "d:o:b:f:s:c:p:a:n:r:l:kh", long_options, &option_index)) != -1) {

		switch (c) {

		// Read in the files
		case 'o':
			output = optarg;
			printf ("saving output to :  %s \n", output.c_str());
			break;
		case 'f':
			fasta = optarg;
			printf ("Parsing fasta file :  %s \n", fasta.c_str());
			break;
		case 'b':
			file = optarg;
			printf ("Using bam file :  %s \n", file.c_str());
			break;
		case 'r':
			ROA_SIZE = atoi(optarg);
			ROA_SIZE = (ROA_SIZE % 2) ? ROA_SIZE+1 : ROA_SIZE;
			if(Read::max_length < ROA_SIZE * 3)
			{
				printf ("ERROR   : max ROA size cannot be set to %d. Set Read::max_length to %d \n", ROA_SIZE, ROA_SIZE*3 );
				ROA_SIZE = Read::max_length / 3;
				printf ("WARNING : max ROA size (b/c of Read::max_length) is %dbps \n", ROA_SIZE);
			}
			printf ("Setting the ROA size to :  %dbps \n", ROA_SIZE);
			break;

		case 'c':
			MIN_COUNT_THRESH = atoi(optarg);
			printf ("Setting the minimum directional count filter threshold to :  %d \n",MIN_COUNT_THRESH);
			break;
		case 'p':
			MIN_FREQ_THRESH = atof(optarg);
			printf ("Setting the minimum directional frequency filter threshold to :  %f \n",MIN_FREQ_THRESH);
			break;
		case 'd':
			DICTIONARY_LEVEL = atoi(optarg);
			printf ("Setting the dictionary verbosity to :  %d \n",DICTIONARY_LEVEL);
			break;

			// Frag making thresholds
		case 'a':
			FRAG_THRESHOLD = atof(optarg);
			printf ("Setting the verified word frequency threshold for fragment assembly to :  %f \n",FRAG_THRESHOLD);
			break;
		case 'n':
			NON_VERIFIED_THRESHOLD = atof(optarg);
			printf ("Setting the non-verified word frequency threshold for fragment assembly to :  %f \n",NON_VERIFIED_THRESHOLD);
			break;

			// SNV Calling parameters
		case 's':
			SNV_VERIFIED_THRESHOLD = atof(optarg);
			printf ("Setting the SNV Verified Call Reason 2 Threshold to :  %f \n",SNV_VERIFIED_THRESHOLD);
			break;
		case 'l':
			SNV_TOTAL_THRESHOLD = atof(optarg);
			printf ("Setting the SNV Total Call Reason 3 Threshold to :  %f \n",SNV_TOTAL_THRESHOLD);
			break;


			// Process Flags
		//case 't':
		//	test();
		//	return EXIT_SUCCESS;
		case 'h':
			usage();
			return EXIT_SUCCESS;
		case 'k':
			DROPID= false;
			printf ("Keeping reads with both inserts and deletes \n");
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

	if( file.empty() || fasta.empty())
	{
		cout << endl << "You must specify a fasta and bam file for DDiMAP to process." << endl << endl;
		usage();
		return EXIT_FAILURE;

	}

	mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


	// ------------------------------------------------------------------------
	// DDiMAP
	// ------------------------------------------------------------------------
	t = clock();
	unique = readFile(file, fasta, ROA_SIZE, DROPID, buildRead);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to read %d | %d reads from BAM file.\n",
			t, ((float)t)/CLOCKS_PER_SEC, unique, total);

	t = clock();
	sequential(MIN_COUNT_THRESH, MIN_FREQ_THRESH, FRAG_THRESHOLD, NON_VERIFIED_THRESHOLD);
	int verified = printFasta(output);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to verify %d | %d.\n",
			t, ((float)t)/CLOCKS_PER_SEC , verified, unique);


	t = clock();
	verified = iterate(buildHistograms);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to build the histogram.\n",
			t, ((float)t)/CLOCKS_PER_SEC);

	//  this code wrote files used for plotting during debugging and is not needed
        //t = clock();
	//printHistograms( output );
	//t = clock() - t;
	//printf ("It took me %lu ticks (%f seconds) to print the histogram.\n",
	//		t, ((float)t)/CLOCKS_PER_SEC);


	t = clock();
	iterate(count);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to count.\n",
			t, ((float)t)/CLOCKS_PER_SEC);


	t = clock();
	callSNVs(SNV_VERIFIED_THRESHOLD, SNV_TOTAL_THRESHOLD, output);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to call SNVs.\n",
			t, ((float)t)/CLOCKS_PER_SEC);

	printDicitonaries(output, DICTIONARY_LEVEL);
//		gnuplot_ctrl    *   h1;
//		h1 = gnuplot_init() ;
//		gnuplot_resetplot(h1) ;
//		gnuplot_cmd(h1, "set xrange [0:1200]");
//		gnuplot_cmd(h1, "set yrange [0.0001:0.9999]");
//		gnuplot_cmd(h1, "set multiplot title \"GC Content for different genes\" layout 3,4");
//		gnuplot_cmd(h1, "plot \"output/Bcl2CG.txt\"     with points ls 1 title \"Bcl2\"");
//		gnuplot_cmd(h1, "plot \"output/Bcl6CG.txt\"     with points ls 1 title \"Bcl6\"");
//		gnuplot_cmd(h1, "plot \"output/CD83CG.txt\"     with points ls 1 title \"CD83\"");
//		gnuplot_cmd(h1, "plot \"output/EmuCG.txt\"      with points ls 1 title \"Emu\"");
//		gnuplot_cmd(h1, "plot \"output/J6-J4CG.txt\"    with points ls 1 title \"j6-j4\"");
//		gnuplot_cmd(h1, "plot \"output/KSCG.txt\"       with points ls 1 title \"KSCG\"");
//		gnuplot_cmd(h1, "plot \"output/Pax5CG.txt\"     with points ls 1 title \"Pax5\"");
//		gnuplot_cmd(h1, "plot \"output/Pim1CG.txt\"     with points ls 1 title \"Pim1\"");
//		gnuplot_cmd(h1, "plot \"output/RhoHCG.txt\"     with points ls 1 title \"RhoH\"");
//		gnuplot_cmd(h1, "plot \"output/Vh1-18CG.txt\" with points ls 1 title \"vh1\"");
//		gnuplot_cmd(h1, "plot \"output/myc-1CG.txt\"    with points ls 1 title \"myc-1\"");
//		gnuplot_cmd(h1, "plot \"output/myc-2CG.txt\" with points ls 1 title \"myc-2\"");
//		gnuplot_cmd(h1, "unset multiplot") ;


	//	test(fasta);

	// ------------------------------------------------------------------------
	// End.
	// ------------------------------------------------------------------------

	return EXIT_SUCCESS;

}
