#ifndef READHELPER_H
#define READHELPER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

using namespace std;

char UINT64ToChar(uint64_t ch, bool upper_case);
string UINT64ToStringCompare(uint64_t s, uint64_t t);
string UINT64ToString(uint64_t s);
uint64_t stringToUINT64(string s);
float CalculateGC(uint64_t seq);
bool hasDash(uint64_t seq);

#endif
