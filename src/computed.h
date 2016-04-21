#ifndef COMPUTED_H

#define COMPUTED_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"


class computed { 
	public:
		computed (int argc, char *argv[]) ;
        int positiontojbin (int chr, int index) ;
		void compute () ;
		void computepairs () ;
		void read_weights (string filename) ;

		data *d;
		data *p;
		genotype *g;


		double maxdis;
		double binsize;
		int nbins;
        bool isgenotype;

		double fcutoff;
		int derived;
		int poly;
		int ascertain;
		int uniform;


		string outfile;
		bool out;
		bool outputpair;
		string outputpairfile;



		unordered_map <string, double> weightmap;
		bool pairs;
		string pairfile;

		bool toflip;
		string flipfile;
		unordered_map <string, string> flipmap;

		string jackpath;
		bool dojack;
        double jackblock;
        int jackminsnps;
        bool jackgenetic;

		string gsfile;
		bool goodsnps;
		bool badsnps;
        bool sample;
		unordered_map<string, string > *goodsnpsmap;
		unordered_map<string, string > *badsnpsmap;


        bool givensnpfilter;
        string filtersnps ; 
        vector<int> cumnumjbins;
};
#endif
