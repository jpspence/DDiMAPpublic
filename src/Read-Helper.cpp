#include "Read-Helper.h"


// ----------------------------------------------------------------------------
// Utitlity Methods
// ----------------------------------------------------------------------------

// Cipher is given by :
uint64_t a    = 0b00000001;
uint64_t c    = 0b00000010;
uint64_t g    = 0b00000011;
uint64_t t 	  = 0b00000100;
uint64_t dash = 0b00000101;

uint64_t charToUINT64(char ch)
{
	switch(ch){
	case 'A':
	case 'a': return a;
	case 'T':
	case 't': return t;
	case 'C':
	case 'c': return c;
	case 'G':
	case 'g': return g;
	case '-': return dash;
	default : return 0;
	}
}

char UINT64ToChar(uint64_t ch, bool upper_case)
{
	if(upper_case){
		if( ch == a) return 'A';
		if( ch == t) return 'T';
		if (ch == c) return 'C';
		if (ch == g) return 'G';
	} else {
		if( ch == a) return 'a';
		if( ch == t) return 't';
		if (ch == c) return 'c';
		if (ch == g) return 'g';
	}
	if (ch == dash) return '-';
	return '\0';
}

string UINT64ToStringCompare(uint64_t s, uint64_t t)
{

	std::stringstream temp;
	while(s!=0){
		temp << UINT64ToChar( s & 0b111 , ( (s & 0b111) == (t & 0b111)) );
		s = s >> 3; t = t >> 3;
	}
	return temp.str();

}

string UINT64ToString(uint64_t s)
{
	return UINT64ToStringCompare(s,s);
}

uint64_t stringToUINT64(string s)
{

	uint64_t temp = 0;
	for ( int i = 0 ; i < s.length()  ;  i++)
		temp += charToUINT64(s[i]) << (3 * i);
	return temp;

}

bool hasDash(uint64_t seq)
{
	while(seq!=0){
		if((seq & 0b111) == dash)
			return true;
		seq = seq >> 3;
	}
	return false;
}

float CalculateGC(uint64_t seq)
{
	float gc = 0.0, total = 0.0;
	while(seq!=0){
		total += 1;
		if((seq & 0b111) == g || (seq & 0b111) ==c )
			gc+=1;
		seq = seq >> 3;
	}
	return gc / total * 100.0;
}
