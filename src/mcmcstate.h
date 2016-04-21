#ifndef MCMC_H

#define MCMC_H

#include "std.h"
#include "printable.h"
#include "data.h"
#include "state.h"


class mcmcstate : public state{ 
	public:

        int iters;
        int burnin;
        int printint;
        //
		// Parameters passed to this program
		data *d;


        mcmcstate ();
        mcmstate (int argc, char *argv[]);

        virtual double propose (mcmcstate &s ) = 0 ;
        void init ();
        void run ();
        void finish ();
}

mcmcstate::mcmcstate () {}


mcmcstate::run () {
    int j  = 0;
    mcmcstate proposed;
    init ();
	for (int i  = 0 ; i  < mcmc_iters; i++) { 
        double q = propose(proposed);
        double u = gsl_ran_flat (rng_r, 0, 1);
        if (u < q ) 
            *this = proposed;
    }
    finish ();
}
#endif
