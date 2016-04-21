#include "mh.h"
#include "io.h"
#include "vectorfn.h"

mh::mh (int argc, char *argv[], posterior *f) {
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
	p->get_int ("debug", io::debug, 0);
	p->get_int ("mcmc_iters", mcmc_iters, 1200);
	p->get_int ("burnin", burnin, 200);
    p->get_int ("thin", thin, 1);
	p->get_int ("print_interval", printint, 1);
	p->get_int ("seed",seed,1,false);

    gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	rng_r = gsl_rng_alloc (rng_T);
	gsl_rng_set (rng_r, seed);

    this->f = f;
    n = f->get_dimension ();
    param.resize (n);
}



void mh::init (vector<double> &initparam) {
    for (int i = 0 ; i < n; i++)
        param[i] = initparam[i];
    logposterior = f->get_posterior (param);
    acceptcount = 0 ;
}

void mh::update ()  {
    vector<double> proposed_param (param);
    double q = f->get_next (param, proposed_param, rng_r);
    double tmppost = f->get_posterior (proposed_param); 
    double acceptprob = 0 ;
    if (std::isnan(tmppost) || std::isinf (tmppost) || std::isnan (q) || std::isinf (q) )
        acceptprob = 0 ;
    else  {
        acceptprob = tmppost + q - logposterior;
        acceptprob = exp(acceptprob);
    }
    cout << tmppost << "\t" << logposterior << "\t" << acceptprob << endl;

    double u = gsl_ran_flat (rng_r,0,1);
    if (u < acceptprob) {
        for (int i = 0 ; i <  n ; i++)
            param[i] = proposed_param[i];
        logposterior = tmppost;
        acceptcount ++;
    }
}


void mh::iterate (vector< vector<double> > &chain) { 
    vector<double> initparam (n,1);
    iterate (initparam, chain);
}

void mh::iterate (vector<double> &initparam, vector< vector<double> > &chain) { 
    init (initparam);
    int total = (mcmc_iters)/thin;
    chain.resize (total, vector<double>(n)); 
    for (int i = 0 ; i  < burnin + mcmc_iters; i++) { 
        update ();
        if (printint > 0 && i % printint == 0 ) { 
			if ( i < burnin) {
                cout << "burnin : " << i << "\t" ; vectorfn::printvector (param) ; cout << "\t" << logposterior << endl;
            } else { 
                cout << "iteration : " << (i-burnin) << "\t" ; vectorfn::printvector (param) ; cout << "\t" << logposterior << endl;
            }

        }
        if (i>=burnin && (i-burnin)%thin == 0 ) {
            int k=(i-burnin)/thin;
            for (int j =  0 ; j  <  n; j++)
                chain[k][j] = param[j] ;
        }
    }
    if (burnin + mcmc_iters > 0)
        acceptrate = (1.*acceptcount)/(burnin + mcmc_iters);
    else
        acceptrate = 0;
    cout << "accept rate = " << acceptrate << endl;
}

