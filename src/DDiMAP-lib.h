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

void iterate(void (*f)(int, int, string, int) );

void print(int gene, int roa, string seq, int count);

void read( BamAlignment ba, int length );

#endif
