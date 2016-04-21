#include "intervalmap.h"
#include "packedgtype.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"
#include "std.h"
#include "functions.h"
#include "fileio.h"
#include "ascertain.h"


ascertain::ascertain (data *d) { 
    astring = "";
    d->get_string ("ascertain", astring, false);
    init (astring);
}

ascertain::ascertain (string s) {
    init (s);
}

void ascertain::init (string s) { 
    vector<string> toks;
    functions::tokenize (s, toks, ";");
    nclauses = toks.size();
    ids.resize (nclauses, vector<string> ()) ;
    nderived.resize (nclauses, vector<int> ()) ;
    ntotal.resize (nclauses, vector<int> ()) ;
    nliterals.resize (nclauses,0);
    for (int i = 0; i < toks.size(); i++) {
        vector<string> t;
        functions::tokenize (toks[i], t, ",");
        nliterals[i]=t.size();
        for (int j = 0; j < t.size(); j++){
            vector<string> t1;
            functions::tokenize (t[j], t1, ":");
            if (t1.size() < 3) {
                cerr << "Ascertain : Badly formatted string " << t[j] << endl;
                exit(1);
            }
            string id = t1[0];
            cout << "id = " << id << endl;
            int a = atoi(t1[1].c_str()); int b = atoi(t1[2].c_str());
            ids[i].push_back (id);
            nderived[i].push_back (a);
            ntotal[i].push_back(b);
        }
    }
}


void ascertain::print ()  {
    for (int i  = 0 ; i < nclauses; i++) { 
        for (int j = 0 ; j < nliterals[i]; j++)  {
            cout << ids[i][j] << "," << nderived[i][j] <<","<<ntotal[i][j] << endl;
        }
    }
}


/*
int main (int argc, char *argv[]) { 
    ascertain a  = ascertain(string(argv[1]));
    a.print ();
}*/
