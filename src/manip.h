#ifndef CONVERTF_H

#define CONVERTF_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"


class manip { 
	public:
		manip (int argc, char *argv[]) ;
		void run ();

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

		// genotype data
		genotype *g;

		bool out;
		string snpoutfile;
		string indivoutfile;
		string genotypeoutfile;
        string inputformat;

		bool uniqpos;
        bool isgenotype;
        bool isx;
        string outputformat;
        string filterpolicy;

        bool extracttrios;
        vector<string> trioind;
        bool trioflag;
};
#endif
