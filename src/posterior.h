#ifndef POSTERIOR_H

#define POSTERIOR_H

#include "std.h"

class posterior { 
	public:
        int proposal_type; 
        vector<double> proposal_params;
        virtual void set_proposal_type (int, vector<double> &) = 0 ;
		virtual double get_posterior (vector<double> &parameter) = 0 ;
        virtual int get_dimension () = 0 ;
        virtual double get_next ( vector<double> & cur, vector<double> & proposed, gsl_rng *) = 0 ;
};

#endif
