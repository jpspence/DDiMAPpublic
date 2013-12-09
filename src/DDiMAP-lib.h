#ifndef DDIMAPLIB_H_GUARD
#define DDIMAPLIB_H_GUARD

#include <map>
#include <utility>
#include "../include/bamtools/src/api/BamAlignment.h"

namespace BamTools {
	struct BamAlignment;
} /* namespace BamTools */

using namespace BamTools;
using namespace std;

struct Read {
	// Each half of a read encodes for up to 32 base pairs
	uint64_t right_sequence_half;
	uint64_t  left_sequence_half;
	unsigned int count;
	unsigned int verification_flags;
};


int iterate ( int (*f)(int, int, string, Read) );
int print (int gene, int roa, string seq, Read read);
int count (int gene, int roa, string seq, Read read);
int verify ( int gene, int roa, string seq, Read read);

void read( BamAlignment ba, int length );

Read buildRead(string word);


#endif
