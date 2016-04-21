#ifndef AFS_H

#define AFS_H

#include "std.h"
#include "printable.h"
#include "data.h"


class afs { 
    public:
        int npops;
        vector<int> chrsperpop;
        vector<string> poplabels;
        unordered_map <string, int> popmap;
        int countsize ; 
        double *counts;
        bool alloced ;
        
        afs (vector<int> &p) ;
        afs (string filename)  ;
        afs ();
        afs (vector<int> &p, vector<string> &pl) ;
        afs* ascertain (vector<string> &pop, vector<int> &dc, vector<int>  &ac) ;
        afs* ascertain (string asc) ;
        afs* permute (vector<string> labels) ;
        afs* permute (string labels);

        ~afs () ;
        double get (vector<int> &index) ;
        void put (vector<int> &index, double v) ;
        void increment (vector<int> &index) ;
        void increment (vector<int> &index, double v) ;
        void increment (int *index, int n ) ;
        void read (string filename ) ;
        void print ()  ;
        void print (string filename);

        double operator () (vector<int> &index ) ;
        double operator () (int *index, int n ) ;


};

#endif

