#ifndef DSTAT_H

#define DSTAT_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"
#include "genomecompute.h"


class f4ratio  : public genomecompute { 
	public:
		f4ratio  (int argc, char *argv[]) ;
        void compute (vector<int> &gtype, int jackindex) ;
        void setparams () ;
        void summarize ();

	    int nstats;
        vector<double> f4rationum ; 
        vector<double> f4ratiodenom ; 
        vector<double> w;
        vector< vector<string> > pops;
        vector< vector< double> > jf4rationum;
        vector< vector< double> > jf4ratiodenom;
        vector< vector< double> > jf4ratiow;

};
#endif
