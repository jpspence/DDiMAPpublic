#ifndef DDIMAPLIB_H

#define DDIMAPLIB_H

#include <map>
#include <utility>
#include <api/BamAlignment.h>
#include <api/BamReader.h>

namespace BamTools {
	struct BamAlignment;
} /* namespace BamTools */

using namespace BamTools;
using namespace std;

struct Read {
	// Each half of a read encodes for up to 32 base pairs
	char sequence[34];
	uint64_t right_sequence_half = 0;
	uint64_t  left_sequence_half = 0;
	unsigned int count;
	unsigned int verification_flags;
};

void readFile(string file, int lenght, Read (*f)(string &, int));
int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) );
int iterate ( int (*f)(int, int, string, Read) );
int print (int gene, int roa, string seq, Read read);
int count (int gene, int roa, string seq, Read read);
int verify ( int gene, int roa, string seq, Read read);
Read convert(string &word, int length);
Read buildRead(string &word, int length);

#endif
