#include "../std.h"

class samplealpha {
	public :
		samplealpha (int n, vector<double> &z, vector<double> &g) {
			this->n = n ; 
			for (int i = 0 ; i  < n; i++){
				(this->z).push_back(z[i]);
				(this->g).push_back(g[i]);
			}
		}

		int n;
		vector<double> z;
		vector<double> g;
		static double  logd ( double x , void * param) {
			double res = 0;
			samplealpha *s = (samplealpha *)param;
			int n = s->n;
			vector<double> &z = s->z;
			vector<double> &g = s->g;

			for (int i = 1 ; i < n ; i++){
				res += log(z[i]) * g[i] - z[i];
			}
			res *=x ;
			for (int i = 1; i < n; i++){ 
				res += (x*log(x)*g[i]);
				res -= gsl_sf_lngamma (x*g[i]);
			}
			return res;
		}

};
