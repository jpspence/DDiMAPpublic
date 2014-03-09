#ifndef DDIMAPLIB_H
#define DDIMAPLIB_H

#include "Read-Helper.h"
#include <map>
#include <utility>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <kseq.h>
#include <api/BamAux.h>
#include <api/SamHeader.h>
#include <api/SamSequence.h>
#include <api/SamSequenceDictionary.h>
#include <kseq.h>
#include <zlib.h>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <regex>
#include <iomanip>
#include <ostream>



namespace BamTools {
struct BamAlignment;
} /* namespace BamTools */

using namespace BamTools;

struct Read {

	// Each half of a read encodes for up to 32 base pairs
	int32_t RefID;
	char sequence[34];
	uint64_t right_sequence_half;
	uint64_t  left_sequence_half;
	unsigned int forward_count;
	unsigned int reverse_count;

	int16_t verification_flags;
	//	0b 0000  0000  0000
	//           ||||  ||||_ 1   Is the read verified on left
	//           ||||  |||__ 2   Is the read verified on right
	//           ||||  ||___ 4   Is this a 50 BP READ (e.g. 50M)
	//           ||||  |____ 8   Does read match reference on left
	//           ||||_		16  Does read match reference on right
	//           |||_ 		32  Does read show up in at least 750ppm
	//           ||_ 		64  Does read show up in at least 1% as a verified frag
	//           |_ 		128 Does read show up in at least 10% non-verified frag

	void set_left_verified()
	{ verification_flags = verification_flags |       0b1; }

	void set_right_verified()
	{ verification_flags = verification_flags |      0b10; }

	void set_left_verified_at_frag()
	{ verification_flags = verification_flags | 0b100000001; }

	void set_right_verified_at_frag()
	{ verification_flags = verification_flags |  0b1000000010; }


	void set_no_indels()
	{ verification_flags = verification_flags |     0b100; }

	void set_matches_ref_on_left()
	{ verification_flags = verification_flags |    0b1000; }

	void set_matches_ref_on_right()
	{ verification_flags = verification_flags |   0b10000; }

	void set_above_ppm_threshold()
	{ verification_flags = verification_flags |  0b100000; }

	void set_above_frag_threshold()
	{ verification_flags = verification_flags | 0b1100000; }

	void set_above_nv_threshold()
	{ verification_flags = verification_flags | 0b11100000; }

	void set_above_nv_threshold_only()
	{ verification_flags = verification_flags | 0b10000000; }

	bool is_left_verified()
	{ return (verification_flags & 0b0000001) ==       0b1;}

	bool is_right_verified()
	{ return (verification_flags & 0b0000010) ==      0b10;}

	bool is_right_left_verified()
	{ return (verification_flags & 0b0000011) ==      0b11;}

	bool is_right_left_verified_at_frag()
	{ return (verification_flags & 0b1100000011) == 0b1100000011;}

	bool matches_ref_on_left()
	{ return (verification_flags & 0b0001000) ==    0b1000;}

	bool matches_ref_on_right()
	{ return (verification_flags & 0b0010000) ==   0b10000;}

	bool matches_reference()
	{ return (verification_flags & 0b0011000) ==   0b11000;}

	bool is_above_ppm()
	{ return (verification_flags & 0b0100000) ==  0b100000;}

	bool is_above_frag()
	{ return (verification_flags & 0b1000000) == 0b1000000;}

	bool is_above_non_verified()
	{ return (verification_flags & 0b10000000) == 0b10000000;}

	unsigned int total_count()
	{ return forward_count + reverse_count;}

};

Read buildRead(string &word, int length);
int readFile(string file, char *fasta, int length, Read (*f)(string &, int));
int iterate ( int (*f)(string, int, string, Read&) );
int printFasta();
void sequential(int threshold, double ppm, double frag, double non_verified);
void callSNVs(double snv_verified_threshold, double snv_total_threshold);
int buildHistograms(string gene, int position, string seq, Read& read);
void printHistograms();
int count (string gene, int position, string seq, Read& read);
void frequency_filter(string gene, int position, int threshold, double ppm, double frag, double non_verified, bool testing, string name, string sequence, int test_position);
void check_verify ( Read r, bool is_right, string gene, int position);
#endif
