#ifndef MERGEIT_H

#define MERGEIT_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"


class mergeit { 
	public:
		mergeit (int argc, char *argv[]) ;
		void merge ()  ;

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

		// genotype data
		genotype *g1;
		snpmap *g2;

		bool out;
		string snpoutfile;
		string indivoutfile;
		string genotypeoutfile;

		bool flipsnps;
		bool badsnps;
		bool addsnps;

		bool keepgen;
};
#endif
