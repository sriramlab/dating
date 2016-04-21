#ifndef GENOMECOMPUTE_H

#define GENOMECOMPUTE_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"


class genomecompute { 
	public:
		genomecompute (int argc, char *argv[]) ;
        int positiontojbin (int chr, int index) ;
		void run ();
        virtual void compute (vector<int> &gtype, int)  = 0 ;
        virtual void setparams () = 0 ;

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

        bool isgenotype;

		string jackpath;
		bool dojack;
        double jackblock;
        int jackminsnps;
        bool outputjack;
        string outputjackstring;
        bool jackgenetic;

        vector<int> cumnumjbins;


        int nind, nchr, njbins;
        
};
#endif
