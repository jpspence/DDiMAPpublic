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
};

int readFile(string file, int length, Read (*f)(string &, int));
int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) );
int iterate ( int (*f)(string, int, string, Read) );
void iterateAndSet( Read reads_array[]);
int print (string gene, int roa, string seq, Read read);
int count (string gene, int roa, string seq, Read read);
int verify ( string gene, int roa, string seq, Read read);
Read convert(string &word, int length);
Read buildRead(string &word, int length);

#endif
