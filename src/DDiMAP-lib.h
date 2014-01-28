#ifndef DDIMAPLIB_H

#define DDIMAPLIB_H

#include <map>
#include <utility>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <iostream>
#include <sstream>
#include <iterator>

namespace BamTools {
struct BamAlignment;
} /* namespace BamTools */

using namespace BamTools;
using namespace std;

struct Read {

	// Each half of a read encodes for up to 32 base pairs
	char sequence[34];
	uint64_t right_sequence_half;
	uint64_t  left_sequence_half;
	unsigned int count;

	unsigned int verification_flags;
	//	0b00000000
	//    ||||||||_	1   Is the read verified on left
	//    |||||||_  2   Is the read verified on right
	//    ||||||_   4   Is this a 50 BP READ (e.g. 50M)
	//    |||||_ 	8   Does read match reference on left
	//    ||||_		16  Does read match reference on right
	//    |||_ 		32  Does read show up in at least 750ppm
	//    ||_ 		64
	//    |_ 		128

	// Function Declarations
	bool matches_ref_on_left_and_right()
	{ return (verification_flags & 0b00011000) == 0b00011000; }

	bool is_right_left_verified()
	{ return (verification_flags & 0b00000011) == 0b00000011; }

	bool is_above_threshold()
	{ return (verification_flags & 0b00100000) == 0b00100000; }

};

int readFile(string file, char *fasta, int length, Read (*f)(string &, int));
int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) );
int iterate ( int (*f)(string, int, string, Read) );
void iterateAndSet( Read reads_array[]);
int print (string gene, int position, string seq, Read read);
int count (string gene, int position, string seq, Read read);
int verify ( string gene, int position, string seq, Read read);
Read convert(string &word, int length);
Read buildRead(string &word, int length);
int buildHistograms(string gene, int position, string seq, Read read);
void printHistograms();
#endif
