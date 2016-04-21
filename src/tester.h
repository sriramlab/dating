#ifndef TESTER_H

#define TESTER_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"


class tester { 
	public:
		tester (int argc, char *argv[]) ;
		void test ();

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

		// genotype data
		genotype *g;

		string snpoutfile;
		string indivoutfile;
		string genotypeoutfile;

};
#endif
