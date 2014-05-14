#ifndef DDIMAPLIB_H
#define DDIMAPLIB_H

#include "Read-Helper.h"
#include <utility>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <iostream>
#include <iterator>
#include <stdio.h>
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
#include <iomanip>
#include <ostream>

namespace BamTools {
struct BamAlignment;
} /* namespace BamTools */

using namespace BamTools;

Read buildRead(string &word, int length);
int readFile(string file, string fasta, int roa_length, bool dropID, Read (*f)(string &, int));
int iterate ( int (*f)(string, int, string, Read&) );
int printFasta(string output);
void sequential(int threshold, double ppm, double frag, double non_verified);
void callSNVs(double snv_verified_threshold, double snv_total_threshold, string output);
int buildHistograms(string gene, int position, string seq, Read& read);
void printHistograms(string output);
void printDicitonaries(string output, int level);
int count (string gene, int position, string seq, Read& read);
void frequency_filter(string gene, int position, int threshold, double ppm, double frag, double non_verified, bool testing, string name, string sequence, int test_position);
void check_verify ( Read r, bool is_right, string gene, int position);
#endif
