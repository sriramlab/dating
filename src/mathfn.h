#ifndef MATHFN_H
#define MATHFN_H


#include "std.h"
#include <sstream>
/*
template <class T>
inline void swap(T &x, T &y) {
	T tmp = x;
	x = y;
	y = x;
}
*/
class mathfn {
	public:
		static int round (double x);
		static double min (double, double);
		static double max (double, double);
		static double lsumexp (double, double);
		static double lminusexp (double, double);
        static double choose (int n, int k ) ;
        static double lnchoose (int n, int k ); 
        static double sigmoid (double);
        static double dsigmoid (double ) ;
};

inline int mathfn::round (double x) { 
	int x1 = (int) floor (x);
	int x2 = (int) ceil (x);
	if (fabs (x-x1) < fabs (x- x2))
		return x1;
	else 
		return x2;
}

inline double mathfn::min(double x, double y) {
	double z = x<y?x:y;
	return z;
}

inline double mathfn::max(double x, double y) {
	double z = x>y?x:y;
	return z;
}

inline double mathfn::lsumexp( double x, double y) {
		double z = max(x,y);
		double w = min(x,y);
		z += log ( 1 + exp(w-z));
		return z;
}

inline double mathfn::lminusexp( double x, double y) {
		if (x < y )  {
				cerr << "Illegal operation in mathfn::lminusexp"<<endl;
				exit(1);
		}
		double z = x + log ( 1 - exp(y-x));
		return z;
}

inline double mathfn::choose (int n, int k ) { 
    if (k>n)
        return 0 ;
    double r = 1;
    k = (k<(n-k))?k:n-k;
    for (int i=1;i<=k;i++) { 
        r *= (n-i+1);
        r/=i;
    }
    return r;
}

inline double mathfn::lnchoose (int n, int k ) { 
    if (k>n)
        return -1.0/0.0 ;
    double r = 0;
    k = (k<(n-k))?k:n-k;
    for (int i=1;i<=k;i++) { 
        r += log(n-i+1);
        r-= log(i);
    }
    return r;
}

inline double mathfn::sigmoid (double x ) { 
    return 1/(1+exp(-x));
}

inline double mathfn::dsigmoid (double x ) { 
    return sigmoid(x)*(1-sigmoid(x));
}
#endif


