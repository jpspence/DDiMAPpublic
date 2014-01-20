#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

using namespace std;

string seq;
map<string, int> patterns;

string BuildSequence(int length)
{

	for(int i = 0; i < length; i++){
		int j = random();
		if(((j & 0b00000011) == 0b00000011))
			seq.append("a");
		else if(((j & 0b00000011) == 0b00000010))
			seq.append("c");
		else if(((j & 0b00000011) == 0b00000001))
			seq.append("t");
		else
			seq.append("g");
	}

	return seq;
}

void CountPatterns(int length)
{
	for(int i = 0; i < seq.length()-length; i++){
		string pattern = seq.substr(i, length);
		patterns[pattern] = (patterns[pattern]) ? patterns[pattern] + 1 : 1;
	}

}

int main (int argc, char **argv)
{
	cout << "hi" <<endl;

	BuildSequence(5000000);
	CountPatterns(22);

	int total = 0;
	map<string, int>::iterator counts;
	for(counts = patterns.begin(); counts != patterns.end(); ++counts){
		if((* counts).second > 1)
				cout<< (* counts).first << " : " << (* counts).second << endl;
		total += (* counts).second;
	}

	cout << "we read " << total <<endl;
}
