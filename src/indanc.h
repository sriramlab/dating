#ifndef DSTAT_H

#define DSTAT_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "genotype.h"
#include "genomecompute.h"


class indanc  : public genomecompute { 
	public:
		indanc  (int argc, char *argv[]) ;
        void compute (vector<int> &gtype, int jackindex) ;
        void setparams () ;
        void summarize ();

	    int nstats;
        vector<double> indancnum;
        vector<double> indancdenom;
        vector<double> w;

        vector< vector<double> > jindancnum;
        vector< vector<double> > jindancdenom;
        vector< vector< double> > jindancw;
        vector<string> pops;


        string freqoutfile;
        bool freqout;
};
#endif
