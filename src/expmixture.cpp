#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "vectorfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "mcmcstate.h"
#include "convertf.h"

ofstream statsfs;


class expmixture: public mcmcstate {
    public:
        int col;
        double l;
        double h;
        double s;
        vector<double> x;
        vector<double> y;
        int n;
        int nmin ;
        int nmax;
        data *p;
 
        bool fitaffine;
        double affine; 

        int ncomponents;
        double sigma2 ; 
        vector<double> intercept;
        vector<double> lambda;
        vector<double> tgf;
        double g[2];


        vector<double> alphav;
        double *alphaprob;
        unsigned int *alphapick;


 
        expmixture ();
        expmixture (int argc, char *argv[]);

        double propose (expmixture &s);
}


expmixture::expmixture (int argc, char *argv[]) {
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
	p = new data (pfile, d->configMap);
	p->print_parameters();
    p->get_int("col",col,2);
    p->get_double("l",l,0.02);
    p->get_double("h",h,1);
    string inpfile;
    p->get_string ("inp",inpfile,true);
    read_rolloff (inpfile);
    string alphafile;
    p->get_string ("alpha", alphafile, true);
    p->get_double ("glb", g[0], 25);
    p->get_double ("gub", g[1], 33);
    string statsfile;
    p->get_string ("output", statsfile, true);
    p->get_double ("tail", tail, 0.05);
    p->get_int ("ncomponents", ncomponents, 1);
    statsfs.open (statsfile.c_str());

    vector<string> *tmp = new vector<string>();
    fileio::read_vector (alphafile, tmp,0);

    alphav = vector<double> (tmp->size());
    alphaprob = new double [tmp->size()];
    alphapick = new unsigned int[tmp->size()];
    for (int i = 0 ; i  < tmp->size();i++){
        alphav[i] = atof((*tmp)[i].c_str());
        alphaprob[i] = 1.0/tmp->size();
        alphapick[i] = 0;
    }

	p->get_int ("debug", io::debug, 0);
	p->get_int ("mcmc_iters", mcmc_iters, 1200);
	p->get_int ("burnin", burnin, 200);
	p->get_int ("print_interval", printint, 1);
	if (burnin > mcmc_iters) {
		cerr << "Number of iterations must be > burnin" << endl;
		exit(1);
	}
	p->get_int ("seed",seed,1,false);

    p->get_boolean ("fitaffine", fitaffine, false, false);

    gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	rng_r = gsl_rng_alloc (rng_T);
	gsl_rng_set (rng_r, seed);

}



double getposterior::logp () {
    double val = 0;
    int np = 0 ;
    double affine = 0  ;
    if ( givenjack ) { 
        val = - 0.5 * pow  ( ( lambda - lambdahat ) , 2) / (lambdase*lambdase) - log  (lambdase);
        return (val);
    } else { 
        if (fitaffine)  {
            int nc = 0 ; 
            for( int i = nmin ; i <= nmax; i++){
                double tmp = y[i];
                for (int j =  0;  j < ncomponents; j++)
                    tmp -= -intercept[j]*exp(-lambda[j]*x[i]/100)
                affine += tmp;
                nc ++;
            }
            affine /= nc;
        } else 
            affine = 0 ;

        for( int i = nmin ; i <= nmax; i++){
            double tmp = y[i];
            for (int j =  0;  j < ncomponents; j++)
                tmp -= -intercept[j]*exp(-lambda[j]*x[i]/100)
            val += tmp*tmp;
            np++;
        }
        val =  -0.5*val/sigma2 - 0.5 * np*log(sigma2) ;
        return (val);
    }
}


double expmixture::propose ( expmixture &prop) { 
    prop = (*this); 
    double q = 0 ;
    for (int j =0 ; j < ncomponents; j++) {
        prop.intercept[j] = exp(gsl_ran_gaussian (rng_r, 0.1) ) * intercept[j];
        prop.lambda[j] = exp(gsl_ran_gaussian (rng_r, 1) ) * lambda[j];
        prop.sigma2  = exp(gsl_ran_gaussian (rng_r, 1)) * sigma2;

        q += (log(prop.intercept[j]) - log(intercept[j]));
        q += (log(prop.lambda[j]) - log(lambda[j]));
        q += (log(prop.sigma2) - log(sigma2));

    }
    q+ = (prop.logp () - logp ());
    return q;

}



void expmixture::init ()  {
}


int main (int argc, char* argv[]) {
    expmixture e = expmixture (argc, argv);
    e.mcmc ();
}
