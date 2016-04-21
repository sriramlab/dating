#ifndef ASCERTAIN_H

#define ASCERTAIN_H
#include "std.h"
#include "snp.h"
#include "printable.h"
#include "data.h"
#include "ind.h"


class ascertain  {
    public:
        int nclauses;
        vector<int> nliterals;
        vector< vector< string> > ids;
        vector< vector< int> > nderived;
        vector< vector< int> > ntotal;
        string astring;

        ascertain (string);
        ascertain (data *d);
        void init (string);
        void print ();
};

#endif
