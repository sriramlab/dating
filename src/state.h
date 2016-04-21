#ifndef STATE_H

#define STATE_H

#include "std.h"
#include "printable.h"
#include "data.h"


class state { 
	public:
        int seed;
		const gsl_rng_type * rng_T;
        gsl_rng * rng_r;

        state (){}
        virtual double logp ()=0;
	
};
#endif
