#ifndef STATSMS_H

#define STATSMS_H

#include "std.h"
#include "printable.h"
#include "data.h"


class statsms1 { 
	public:
		statsms1 (string);
		void readms (string filename);  
		void readmshot (int reps);
;
        int segsub( int nsub, int segsites, char **list );
        int frequency( char allele,int site,int nsam,  char **list);
        double thetah( int nsam, int segsites, char **list);
        double hfay( int nsam, int segsites, char **list);
        double nucdiv( int nsam, int segsites, char **list);

		string paramfile;
		string msfile;
		data *param;

		data *d;
		int length;
		double rho;
		int reps;
		vector<int> samples;
		int nsamples;
        vector<string> samplenames;


		double maxdis;
		double binsize;
		int nbins;

		double fcutoff;
		int derived;
		int poly;
		int uniform;

		bool out;
		string outfile;
		bool fast;
		bool asn;
		bool kg;

		int seed;
		int usegenetic;
		int start;
		bool mshot;
		int corr;


        int ooa;
        int afr;
        int den;
        int nea;
        int chimp;

        unordered_map <string,int> samplemap;

		const gsl_rng_type * T;
	        gsl_rng * r;
		double bgrho;
		string snpfile;
		bool printsnplist;
		double alpha;
        bool hetascertain;
        string ascertain;
        int replicates;
        int chrsperrep;
};
#endif
