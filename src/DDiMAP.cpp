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
#include "DDiMAP-test.h"
#include "../include/gnuplot_i/gnuplot_i.h"

#define SLEEP_LGTH  1

// Default file.
string file  = "/Users/androwis/Dropbox/documents/DataExchangeUR/128_Gen7_CUSHAW2_VhJ_sorted.bam";
string fasta = "/Users/androwis/Dropbox/documents/DataExchangeUR/128_Gen7_CUSHAW2_VhJ.fa";
string output = "./output/";

// Default values
int ROA_SIZE 		   = 34;
int  VERIFY_THRESHOLD  =  2;
bool DROPID = true;

// Frag making thresholds
double PPM = 0.00075;
double FRAG_THRESHOLD = .01;
double NON_VERIFIED_THRESHOLD = .1;

// SNV calling thresholds
double SNV_VERIFIED_THRESHOLD = .003;
double SNV_TOTAL_THRESHOLD = .1;


// Output Parameters
int DICTIONARY_LEVEL = 0;

void usage()
{
	cout << "usage : DDiMAP [-f <fasta> -b <bam> <args>] [--help]" << endl;

	cout << endl << "Basic Parameters:" << endl;
	cout << "   --bam              | -b   This specifies the path to the bam file" << endl;
	cout << "   --fasta            | -f   This specifies the path to the fasta file" << endl;
	cout << "   --keepID           | -k   Keep reads that have both an insert and delete in CIGAR string" << endl;
	cout << "   --verify-threshold | -v   Minimum number of reads to see in each direction (default : " << VERIFY_THRESHOLD<< ")" << endl;
	cout << "   --roa-size         | -r   Number of base pairs for a Region of Analysis    (default : " << ROA_SIZE << ")" << endl;

	cout << endl << "Frag Making Parameters" << endl;
	cout << "   --ppm              | -p   Minimum level of reads to consider for DDiMAP    (default : "<< PPM << ")" << endl;
	cout << "   --frag-threshold   | -t   Minimum verified coverage required to be considered for frags (default : "<< FRAG_THRESHOLD <<")" << endl;
	cout << "   --nv-threshold     | -n   Minimum non-verified coverage required to be considered for frags (default : "<< NON_VERIFIED_THRESHOLD <<")" << endl;

	cout << endl << "SNV Calling Parameters" << endl;
	cout << "   --snv-verified     | -s   Minimum level of nucleotide variation in verified words to call an SNV (default : " << SNV_VERIFIED_THRESHOLD<<")" << endl;
	cout << "   --snv-total        | -l   Minimum level of nucleotide variation in total to call an SNV (default : " << SNV_TOTAL_THRESHOLD<< ")" << endl;

	cout << endl << "Output Parameters" << endl;
	cout << "   --output           | -o   Directory to store output (default : "<< output <<")" << endl;
	cout << "   --dictionary-level | -d   Dictionary verbosity : 0 = fwd/rev counts | 1 = in/del data | 2 = frag mappings (default : "<<DICTIONARY_LEVEL<<")" << endl;

	cout <<endl;
	cout << "Future Parameters (works in progress):"<<endl;
	cout << "   --length-of-snv-ref| -l   Number of base pairs you'd like to see in SNV" << endl;

}
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
			{"verify-threshold", 	0, 0, 'v'},
			{"roa-size", 	0, 0, 'r'},
			{"output", 0,0, 'o'},

			{"ppm", 	0, 0, 'p'},
			{"frag-threshold", 	0, 0, 't'},
			{"nv-threshold", 	0, 0, 'n'},

			{"snv-threshold", 	0, 0, 's'},
			{"snv-total", 	0, 0, 'l'},

			{"help", 	0, 0, 'h'},
			{NULL, 		0, NULL, 0}
	};

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "d:o:b:f:t:s:v:p:n:r:l:kh", long_options, &option_index)) != -1) {

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

		case 'v':
			VERIFY_THRESHOLD = atoi(optarg);
			printf ("Setting the Verify Threshold to :  %d \n",VERIFY_THRESHOLD);
			break;
		case 'd':
			DICTIONARY_LEVEL = atoi(optarg);
			printf ("Setting the Dictionary Verbosity to :  %d \n",DICTIONARY_LEVEL);
			break;

			// Frag making thresholds
		case 'p':
			PPM = atof(optarg);
			printf ("Setting the PPM threshold to :  %f \n",PPM);
			break;
		case 't':
			FRAG_THRESHOLD = atof(optarg);
			printf ("Setting the frag threshold to :  %f \n",FRAG_THRESHOLD);
			break;
		case 'n':
			NON_VERIFIED_THRESHOLD = atof(optarg);
			printf ("Setting the non-verified threshold to :  %f \n",NON_VERIFIED_THRESHOLD);
			break;

			// SNV Calling parameters
		case 's':
			SNV_VERIFIED_THRESHOLD = atof(optarg);
			printf ("Setting the SNV Verified Threshold to :  %f \n",SNV_VERIFIED_THRESHOLD);
			break;
		case 'l':
			SNV_TOTAL_THRESHOLD = atof(optarg);
			printf ("Setting the SNV Total Threshold to :  %f \n",SNV_TOTAL_THRESHOLD);
			break;


			// Process Flags
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
	sequential(VERIFY_THRESHOLD, PPM, FRAG_THRESHOLD, NON_VERIFIED_THRESHOLD);
	int verified = printFasta(output);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to verify %d | %d.\n",
			t, ((float)t)/CLOCKS_PER_SEC , verified, unique);


	t = clock();
	verified = iterate(buildHistograms);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to build the histogram.\n",
			t, ((float)t)/CLOCKS_PER_SEC);

	t = clock();
	printHistograms( output );
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to print the histogram.\n",
			t, ((float)t)/CLOCKS_PER_SEC);


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
	//	gnuplot_ctrl    *   h1;
	//	h1 = gnuplot_init() ;
	//	gnuplot_resetplot(h1) ;
	//	gnuplot_cmd(h1, "set xrange [0:1200]");
	//	gnuplot_cmd(h1, "set yrange [0.0001:0.9999]");
	//	gnuplot_cmd(h1, "set multiplot title \"GC Content for different genes\" layout 3,4");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Bcl2CG.txt\"     with points ls 1 title \"Bcl2\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Bcl6CG.txt\"     with points ls 1 title \"Bcl6\"");
	//	gnuplot_cmd(h1, "plot \"/Users/pa/Desktop/CD83CG.txt\"     with points ls 1 title \"CD83\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/EmuCG.txt\"      with points ls 1 title \"Emu\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/J6-J4CG.txt\"    with points ls 1 title \"j6-j4\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/KSCG.txt\"       with points ls 1 title \"KSCG\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Pax5CG.txt\"     with points ls 1 title \"Pax5\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Pim1CG.txt\"     with points ls 1 title \"Pim1\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/RhoHCG.txt\"     with points ls 1 title \"RhoH\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/Vh1-18CG.txt\" with points ls 1 title \"vh1\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/myc-1CG.txt\"    with points ls 1 title \"myc-1\"");
	//	gnuplot_cmd(h1, "plot \"/Users/androwis/Desktop/myc-2CG.txt\" with points ls 1 title \"myc-2\"");
	//	gnuplot_cmd(h1, "unset multiplot") ;
	//

	//	test(fasta);

	// ------------------------------------------------------------------------
	// End.
	// ------------------------------------------------------------------------

	return EXIT_SUCCESS;

}
