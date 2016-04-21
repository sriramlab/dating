#ifndef DSTAT_H

#define DSTAT_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"
#include "genomecompute.h"


class dstat  : public genomecompute { 
	public:
		dstat  (int argc, char *argv[]) ;
        void compute (vector<int> &gtype, int jackindex) ;
        void setparams () ;
        void summarize ();

	    int nstats;
        vector<double> baba;
        vector<double> abba;
        vector<double> dstatnum ; 
        vector<double> dstatdenom ; 
        vector<double> dw;
        vector< vector<string> > pops;
        vector< vector< double> > jdstatnum;
        vector< vector< double> > jdstatdenom;
        vector< vector< double> > jdstatw;

};
#endif
