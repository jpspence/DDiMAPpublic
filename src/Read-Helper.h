#ifndef READHELPER_H
#define READHELPER_H

#include <cstdint>
#include <map>
#include <sstream>
#include <bitset>


using namespace std;

#include "Read-Helper.h"

struct Read {

	map<int, int> frag_counts;
	map<int32_t, uint64_t> RefID;

	// 0 : no indels
	// 1 : inserts
	// 2 : deletes
	// 3 : in/dels
	uint64_t cigar_counts[4] = {0};

	// This is the number of bits we're using.
	// ROA_SIZE = max_length / 3
	// e.g.  if you want an ROA =  34 ==> max_length = 3
	static int const max_length = 300;
	static int const half_length = max_length/2;
	char sequence[max_length];
	std::bitset<max_length/2> right_sequence_half;
	std::bitset<max_length/2> left_sequence_half;
	unsigned int forward_count = 0;
	unsigned int reverse_count = 0;

	int16_t verification_flags;
	//	0b 0000  0000  0000
	//           ||||  ||||_ 1  Is the read verified on left
	//           ||||  |||__ 2  Is the read verified on right
	//           ||||  ||___ 4  Does read match reference on left
	//           ||||  |____ 8  Does read match reference on right
	//           ||||
	//           ||||_____  16
	//           |||______  32  Does read show up in at least 750ppm
	//           ||_______  64  Does read show up in at least 1% as a verified frag
	//           |________ 128 Does read show up in at least 10% non-verified frag

	void set_left_verified()
	{ verification_flags = verification_flags |       0b1; }

	void set_right_verified()
	{ verification_flags = verification_flags |      0b10; }

	void set_left_verified_at_frag()
	{ verification_flags = verification_flags |  0b100000001; }

	void set_right_verified_at_frag()
	{ verification_flags = verification_flags |  0b1000000010; }

	void set_data(bool hasDeletions, bool hasInsertions, bool isReverseStrand, int32_t refID)
		{
			RefID[refID]++;

			if(isReverseStrand)
				reverse_count++;
			else
				forward_count++;

			if(!hasDeletions && !hasInsertions)
				cigar_counts[0]++;
			else if(hasInsertions && hasDeletions)
				cigar_counts[3]++;
			else if(hasInsertions)
				cigar_counts[1]++;
			else if(hasDeletions)
				cigar_counts[2]++;

		}

	void set_matches_ref_on_left()
	{ verification_flags = verification_flags |    0b100; }

	void set_matches_ref_on_right()
	{ verification_flags = verification_flags |   0b1000; }

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

	bool is_right_verified_at_frag()
	{ return (verification_flags & 0b1100000010) == 0b1100000010;}

	bool is_left_verified_at_frag()
	{ return (verification_flags & 0b1100000001) == 0b1100000001;}

	bool matches_ref_on_left()
	{ return (verification_flags & 0b0001000) ==    0b100;}

	bool matches_ref_on_right()
	{ return (verification_flags & 0b0010000) ==   0b1000;}

	bool matches_reference()
	{ return (verification_flags & 0b0011000) ==   0b1100;}

	bool is_above_ppm()
	{ return (verification_flags & 0b0100000) ==  0b100000;}

	bool is_above_frag()
	{ return (verification_flags & 0b1000000) == 0b1000000;}

	bool is_above_non_verified()
	{ return (verification_flags & 0b10000000) == 0b10000000;}

	unsigned int total_count()
	{ return forward_count + reverse_count;}

};

#endif
