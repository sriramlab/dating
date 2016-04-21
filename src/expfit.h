#ifndef EXPFIT_H

#define EXPFIT_H
#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "vectorfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "convertf.h"

class expfit  { 
    public:
        double findmle(double &mleintercept, double &mlelambda, double &mlesigma2) ;
        void read_input (string filename) ;
        void set_data ( vector<double> &tx, vector<double> &ty) ;
        expfit (int col  = 1, double l = 0, double h = 1, double scale = 1); 

        int col;
        double l;
        double h;
        double scale ;
        int nmin ;
        int nmax;
        vector<double> x;
        vector<double> y;

};


#endif

