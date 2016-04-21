#include "ram.h"
#include "afs.h"

int main (int argc, char *argv[]) { 
    if (argc>1){
        afs a (argv[1]);

//        vector<string> ascpops(1); ascpops [0] =  "y";;
        //    vector<string> ascpops(2); ascpops [0] =  "d"; ascpops[1] = "n";
//        vector<int> ad(1); ad[0] =  1; 
//        vector<int> at(1); at[0] =  2;
//        afs *b = a.ascertain  ( ascpops, ad , at);
//        afs *b = a.ascertain ("ascertain:  y :: 1:2");
//        b->print ();

            afs *b = a.permute ("sampsizes: y: 3; n:1; d:1");
            b->print ();
    }
}
