#ifndef DDIMAPLIB_H

#define DDIMAPLIB_H

#include <map>
#include <utility>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>

namespace BamTools {
struct BamAlignment;
} /* namespace BamTools */

using namespace BamTools;
using namespace std;

struct Read {

	// Each half of a read encodes for up to 32 base pairs
	int32_t RefID;
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

	void set_left_verified()
	{ verification_flags = verification_flags | 0b0000001; }

	void set_right_verified()
	{ verification_flags = verification_flags | 0b0000010; }

	void set_no_indels()
	{ verification_flags = verification_flags | 0b0000100; }

	void set_matches_ref_on_left()
	{ verification_flags = verification_flags | 0b0001000; }

	void set_matches_ref_on_right()
	{ verification_flags = verification_flags | 0b0010000; }

	void set_above_ppm_threshold()
	{ verification_flags = verification_flags | 0b0100000; }

	bool is_right_left_verified()
	{ return (verification_flags & 0b0000011) == 0b0000010;}

	bool matches_reference()
	{ return (verification_flags & 0b0011000) == 0b0011000;}

	bool matches_ref_on_left()
	{ return (verification_flags & 0b0001000) == 0b0001000;}

	bool matches_ref_on_right()
	{ return (verification_flags & 0b0010000) == 0b0010000;}

	bool is_above_ppm()
	{ return (verification_flags & 0b0100000) == 0b0100000;}

	bool right_half_matches_track_left_offset(uint64_t other_track)
	{
		return (right_sequence_half & 0b111111111111111111111111) == (other_track >> 27) ;
	}

	bool right_half_matches_track_left_offset2(uint64_t other_track)
	{
		return (right_sequence_half & 0b111111111111111111111111111) == (other_track >> 24) ;
	}

	bool right_half_matches_track_right_offset(uint64_t other_track){
		return (right_sequence_half >> 24) == (other_track & 0b111111111111111111111111111 ) ;
	}

	bool right_half_matches_track_right_offset_2(uint64_t other_track){
		return (right_sequence_half >> 27) == (other_track & 0b111111111111111111111111 ) ;
	}


	bool matches_track_left_offset(uint64_t other_track)
	{
		return (left_sequence_half & 0b111111111111111111111111111) == (other_track >> 24) ;
	}

	bool matches_track_right_offset(uint64_t other_track){
		return (left_sequence_half >> 27) == (other_track & 0b111111111111111111111111 ) ;
	}


	bool matches_track_right_offset_2(uint64_t other_track){
		return (left_sequence_half >> 24) == (other_track & 0b111111111111111111111111111) ;
	}

};

int readFile(string file, char *fasta, int length, Read (*f)(string &, int));
int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) );
int iterate ( int (*f)(string, int, string, Read&) );
void iterateAndSet( Read reads_array[]);
int print (string gene, int position, string seq, Read& read);
int count (string gene, int position, string seq, Read& read);
int verify ( string gene, int position, string seq, Read& read);
int callSNVs( string gene, int position, string seq, Read& read);
Read convert(string &word, int length);
Read buildRead(string &word, int length);
int buildHistograms(string gene, int position, string seq, Read& read);
void printHistograms();
#endif
