#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "genotype.h"


class couple { 
    public:
      genotype *geno;
      data *d;
      data *p;
      genotype *g;
      bool pack;
      int nind;
      int nsnps;
      int nchr;
    bool isgenotype;
    bool givenindivfilterfile;
        int seed;
		const gsl_rng_type * rng_T;
        gsl_rng * rng_r;

        double distance (int k , int l ) ;
        couple (int argc, char *argv[] ) ;
        void process () ; 
};




couple::couple (int argc, char *argv[] ) { 
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
    p->get_int ("debug", io::debug, 0);

	string geno1, snp1, ind1;
    string inputformat = "eigenstrat";
	p->get_string ("genotypename",geno1,true);
	p->get_string ("snpname",snp1,true);
	p->get_string ("indivname",ind1,true);
	p->get_boolean ("isgenotype", isgenotype, false);
    p->get_boolean ("pack", pack, true, false);
    string indivfilterfile = "";
    givenindivfilterfile = p->get_string ("indivfilterfile", indivfilterfile, false);

    g = new genotype (snp1,ind1, geno1, "eigenstrat", isgenotype, indivfilterfile, "", false, pack);

	p->get_int ("seed",seed,1,false);

    nind = g->nind;
    nsnps = g->nsnp;
    nchr = g->nchr;
    gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	rng_r = gsl_rng_alloc (rng_T);
	gsl_rng_set (rng_r, seed);


}

double couple::distance (int k , int l ) { 
    double d = 0;
	for (int i = 0 ; i < nchr; i++){
        vector<snp> &s = g->snps[g->chrs[i]];
		for (int j  = 0 ; j < s.size();j++) {
            int a = (*g)(i,j,k); int b =(*g)(i,j,l);
            d += (a!=b);
        }
    }
    return d;
}

void couple::process ()  { 
    int nrand = 100;
    double *dobs =new double[nind/2];
    double *drandom = new double [nrand];

    for (int i = 0  ;  i < g->nind; i+=2) { 
        dobs[i/2] = distance (i,i+1);
    }
    for (int i  = 0 ; i < nrand; i++){  
        int j =  gsl_rng_uniform_int (rng_r, nind);
        int k = j;
        do {
            k =  gsl_rng_uniform_int (rng_r, nind);
        } while (k==j);
        drandom [i] = distance (j,k);
    }
    double mobs = gsl_stats_mean (dobs, 1 ,nind/2); double sobs = sqrt (gsl_stats_variance (dobs,1,nind/2));
    double mrandom = gsl_stats_mean (drandom, 1 ,nrand); double srandom = sqrt (gsl_stats_variance (drandom,1,nrand));

    cout << mobs << "\t" << sobs << endl;
    cout << mrandom << "\t" << srandom << endl;
}



int main (int argc, char *argv[]) {
    couple c(argc,argv);
    c.process ();
}

