#ifndef MH_H

#define MH_H
#include "posterior.h"
#include "data.h"

class mh {  

    public:
        data *d, *p;
        int seed;
		const gsl_rng_type * rng_T;
        gsl_rng * rng_r;
        int mcmc_iters;
        int burnin;
        int thin;
        int printint;
        double acceptrate ;
        int acceptcount  ;


        int n ;
        posterior *f;
        vector<double> param;
        double logposterior;

        mh (int argc, char *argv[], posterior *f) ;
        void init (vector<double> &) ;
        void update ()  ;
        void iterate  (vector< vector<double> > &) ; 
        void iterate (vector<double> &initparam, vector< vector<double> > &chain) ;

};
#endif
