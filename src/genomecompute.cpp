#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "genomecompute.h"


genomecompute::genomecompute (int argc, char *argv[]) {
	static const char *optString = "vh";

	static const struct option longOpts[] = {
		{ "parameter", required_argument, NULL, 'p' },
		{ "verbose", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, no_argument, NULL, 0 }
	};
	d = new data (argc, argv, optString, longOpts);
	string pfile;
	d->get_string ("parameter",pfile,true);
	p = new data (pfile);
	p->print_parameters();

    g = new genotype (p);
//    dojack = p->get_string ("jackknife", jackpath, false);
    dojack = true;
    p->get_boolean ("jackgenetic", jackgenetic, false, false);
    if (jackgenetic) {
        p->get_double ("jackblock", jackblock, 5, false);
        jackblock /=100;
    } else {
        p->get_double ("jackblock", jackblock, 10, false);
        jackblock = 1e6 * jackblock;
    }
    jackminsnps = p->get_int ("jackminsnps", jackminsnps, 100, false);
    outputjackstring = "";
    outputjack = p->get_string ("outputjack", outputjackstring, false);
    p->get_int ("debug" , io::debug, 0 );

	
}

int genomecompute::positiontojbin (int chr, int index) { 
    if (jackblock <= 0){
        cerr << "Invalid jackknife block size" << endl;
        exit(1);
    }
    vector<snp> &tmpsnps =  g->snps[g->chrs[chr]];
    int i ;
    double pos;
    if (jackgenetic ) {
        pos = tmpsnps[index].genpos;
        i = (pos - g->mingpos[g->chrs[chr]])/jackblock;
    } else { 
        pos = tmpsnps[index].physpos;
        i = (pos - g->minppos[g->chrs[chr]])/jackblock;
    }
    i = i + cumnumjbins[chr];
    return i;    
}

void genomecompute::run ()  {
	nind = g->nind;
	nchr = g->nchr;
    njbins = nchr;
    cout << "# " <<nind << "\t" << nchr << "\t" << g->nsnp << endl;

    cumnumjbins.resize (nchr, 0);
    if (jackblock > 0 ){
        njbins =  0;
        for (int i = 0 ; i < nchr ; i++){
            int tmp = ceil( (g->maxppos[g->chrs[i]] - g->minppos[g->chrs[i]])/jackblock) ;
            if (jackgenetic) 
                tmp = ceil( (g->maxgpos[g->chrs[i]] - g->mingpos[g->chrs[i]])/jackblock) ;
            cumnumjbins[i] = njbins;
            njbins += tmp;
        }
    } else {
        njbins = nchr;
    }
    cout << "# Number of jacknknife blocks = " << njbins << endl;
    setparams ();

    vector <int> gtype (g->nind, 0);
    for (int i = 0 ; i < nchr; i++){
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];
		for (int j  = 0 ; j < tmpsnps.size();j++) {
			snp &s = tmpsnps[j];
            g->get_geno (i, j, gtype);
            int jj = i;
            if (dojack) {
                if (jackblock >0 ) {
                    jj = positiontojbin (i, j);
                } 
            }
            if (io::debug >= 2)
                cout << s.id << "\t" << jj << endl;

            compute (gtype, jj);
        }
    }
}

