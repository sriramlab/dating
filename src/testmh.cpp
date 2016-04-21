#include "posterior.h"
#include "mh.h"
#include "vectorfn.h"

class gammaposterior:public posterior {
    public:
        gammaposterior () { 
            this->proposal_type = 0 ;
            this->proposal_params["sd"] = 1;
        }

        gammaposterior (int proposal_type, unordered_map <string, double> &params) { 
            this->proposal_type = proposal_type;
            this->proposal_params = params;
        }

        void set_proposal_type (int proposal_type, unordered_map <string, double> &params) {
            this->proposal_type = proposal_type;
            this->proposal_params = params;
        }

        double get_posterior ( vector<double> &parameter){ 
            return gsl_ran_gamma_pdf (parameter[0], 1.7, 1/4.4);
        }

        int get_dimension () {return 1;}

        double get_next (vector<double> &cur, vector<double> &proposed, gsl_rng *r) {
            switch (proposal_type) { 
                case 0:
                    // Symmetric gaussian
                    double sd = (proposal_params.find("sd")!=proposal_params.end())?proposal_params["sd"]:1;
                    for (int i = 0 ; i < cur.size(); i++) {
                        double p = gsl_ran_gaussian (r, 1)  + cur[i] *  sd;
                        proposed[i] = p;
                    }
                    return 1;
                break;
            }
        }
    
};

int main (int argc, char *argv[]) { 
    gammaposterior* gp = new gammaposterior (); 
    unordered_map <string, double> params ; params ["sd"] = 0.1;
    gp->set_proposal_type (0, params);
    mh m (argc, argv, gp);
    vector< vector<double > > chain;
    m.iterate (chain);

    for (int i = 0 ; i  < chain.size(); i++) { 
        vectorfn::printvector (chain[i]);
    }
}

