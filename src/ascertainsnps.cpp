#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "ascertain.h"
#include "genotype.h"

class ascertainsnps : public ascertain {
    public:
		// Parameters specific to the data
		data *d;

		// genotype data
		genotype *g;

        bool nomissing;
        bool onlypolymorphic;

        int seed;
		const gsl_rng_type * rng_T;
        gsl_rng * rng_r;

        ascertainsnps (data *d) ;
        void run ()  ;
};


ascertainsnps::ascertainsnps (data *d) :ascertain (d) { 
    this->d = d;
    g = new genotype (d);    

	d->get_int ("seed",seed,1,false);
    d->get_boolean ("nomissing", nomissing, false);
    d->get_boolean ("onlypolymorphic", onlypolymorphic, false);
    d->get_int ("debug" , io::debug, 0 );
    gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	rng_r = gsl_rng_alloc (rng_T);
	gsl_rng_set (rng_r, seed);

}


void ascertainsnps::run ()  {
    string snpoutfile = "";
	d->get_string("snpoutname|snpoutfilename",snpoutfile,true);
	ofstream snpfs (snpoutfile.c_str());
	int nind = g->nind;
	int nchr = g->nchr;
    cout << "debug = " << io::debug << endl;

   	for (int i = 0 ; i < nchr; i++){
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];
		for (int j  = 0 ; j < tmpsnps.size();j++) {
			snp &s = tmpsnps[j];

            int vc, rc, mc;
            g->get_freq ( i, j, vc, rc, mc,"");
            if (nomissing && mc > 0)
                continue;

            if (onlypolymorphic && (vc==0 || rc==0))
                continue;

            if (astring.compare("")==0) { 
                snpfs << s.id << endl;
            } else {
                bool flag = false;
                if (io::debug>=1 ) {
                    cout <<s.id; 
                }
                for (int k  = 0 ; k < nclauses; k++) { 
                    int l = 0;
                    for ( l = 0 ; l < nliterals[k]; l++)  {
                        int vc, rc, mc;
                        int nd = nderived[k][l];int nt = ntotal[k][l];
                        g->get_freq ( i, j, vc, rc, mc, ids[k][l]);
                        if (io::debug>=1 ) {
                            cout << "\t" << vc << "," << rc << "," << mc;
                        }
                        if (nomissing && mc > 0)
                            break;
                        if (vc + rc == nt){
                            if (vc != nd)
                                break;
                        } else if (vc + rc < nt) {
                            break;
                        } else {
                            int a = gsl_ran_hypergeometric (rng_r, vc, rc, nt);
                            if ( a != nd)
                                break;
                        }
                    }
                    if (l==nliterals[k]) {
                        flag = true; break;
                    }
                }
                if (io::debug>=1){ 
                    cout << endl;
                }

                if (flag) {
                    snpfs << s.id << endl;
                }
            }

        }
    }
    snpfs.close ();
	
}



int main (int argc, char* argv[]) {
    static const char *optString = "vh";

    static const struct option longOpts[] = {
        { "parameter", required_argument, NULL, 'p' },
        { "verbose", no_argument, NULL, 'v' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };
    data *d = new data (argc, argv, optString, longOpts);
    string pfile;
    d->get_string ("parameter",pfile,true);
    data *p = new data (pfile, d->configMap);
    p->print_parameters();


	ascertainsnps a = ascertainsnps(p);
    a.print ();
	a.run ();
}
