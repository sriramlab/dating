#ifndef DSTAT_H

#define DSTAT_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"
#include "genomecompute.h"


class f4  : public genomecompute { 
	public:
		f4  (int argc, char *argv[]) ;
        void compute (vector<int> &gtype, int jackindex) ;
        void setparams () ;
        void summarize ();

	    int nstats;
        bool normalize;
        string normpop;

        vector<double> baba;
        vector<double> abba;
        vector<double> f4num ; 
        vector<double> f4denom ; 
        vector<double> dw;
        vector< vector<string> > pops;
        vector< vector< double> > jf4num;
        vector< vector< double> > jf4denom;
        vector< vector< double> > jf4w;

};
#endif
