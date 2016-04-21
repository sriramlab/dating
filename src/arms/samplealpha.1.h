#include "../std.h"

class samplealpha {
	public :
		samplealpha (int n, double *z, double *g) {
			this->n = n ; 
			this->z = new double [n]; this->g = new double [n];
			for (int i = 0 ; i  < n; i++){
				this->z[i] = z[i];
				this->g[i] = g[i];
			}
		}

		int n;
		double* z;
		double* g;
		static double  logd ( double x , void * param) {
			double res = 0;
			samplealpha *s = (samplealpha *)param;
			int n = s->n;
			double *z = s->z;
			double *g = s->g;

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
