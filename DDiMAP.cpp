//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "bamtools/src/api/BamReader.h"
#include "bamtools/src/api/BamWriter.h"

using namespace BamTools;

int main(void) {
	BamReader *br = new BamReader();
	BamAlignment ba;
	br->Open("data/128test_Gen1_example_sorted.bam");

	while(br->GetNextAlignment(ba))
		if(ba.Position > 0)
			printf("%40s %5d %s \n",ba.Name.c_str(), ba.Position, ba.AlignedBases.c_str());

	puts("Boom.");
	return EXIT_SUCCESS;

}
