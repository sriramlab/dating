#ifndef COMPUTEDMS_H

#define COMPUTEDMS_H

#include "std.h"
#include "printable.h"
#include "data.h"


class computedms { 
	public:
		computedms (string);
		void readms (string filename);  
		void readmshot (int reps);

		string paramfile;
		string msfile;
		data *param;

		data *d;
		int length;
		double rho;
		int reps;
		vector<int> samples;
		int nsamples;


		double maxdis;
		double binsize;
		int nbins;

		double fcutoff;
		int derived;
		int poly;
		int ascertain;
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



		const gsl_rng_type * T;
	        gsl_rng * r;
		double bgrho;
		string snpfile;
		bool printsnplist;
		double alpha;
};
#endif
