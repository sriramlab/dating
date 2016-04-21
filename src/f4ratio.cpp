#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "f4ratio.h"

f4ratio::f4ratio (int argc, char *argv[]) 
    : genomecompute (argc, argv)  {

}


void f4ratio::setparams () { 
    string poplist = "";
    p->get_string ("popfilename", poplist, true);

	ifstream inp (poplist.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< poplist <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
        vector<string> toks;
        functions::tokenize (line.c_str(), toks, " \t");
        pops.push_back (toks);
    }
    inp.close ();
    nstats = pops.size();
    abba.resize (nstats);
    baba.resize (nstats);

    f4rationum.resize (nstats);
    f4ratiodenom.resize (nstats);
    w.resize (nstats);
    jf4rationum.resize ( nstats, vector<double> (njbins) ) ;
    jf4ratiodenom.resize ( nstats, vector<double> (njbins) ) ;
    jf4ratiow.resize (nstats, vector<double> (njbins) );
    for (int i  = 0  ; i < nstats; i ++ ) {
        vectorfn::printvector(pops[i]); cout << endl;
    }

}

void f4ratio::compute (vector<int> &gtype, int jackindex) {
    for (int i  = 0  ; i < nstats; i++ ) { 
        vector<string> &p = pops[i];
        double a = g->get_freq (gtype,  p[0]); 
        double o = g->get_freq (gtype,  p[1]); 
        double x = g->get_freq (gtype,  p[2]); 
        double c = g->get_freq (gtype,  p[3]); 
        double b = g->get_freq (gtype,  p[5]);
        if (a<0 || o<0 || x<0 || c<0 || b <0) {
           continue;        
        }

        double w = (1-c)*d;
        double num = (a-o) * (x-c);
        double denom = (a-o) *(b-c);
        f4rationum[i] += num;
        f4ratiodenom[i] += denom;
        w[i] += 1;
        jf4rationum[i][jackindex] += num;
        jf4ratiodenom[i][jackindex] += denom;
        jf4ratiow[i][jackindex] += 1;
    }

}


void f4ratio::summarize () {
    cout << "#Populations";
    cout << "\tF4 ratio\tBias-corrected F4 ratio\tSE\tZ-score";
    cout << endl;
    for (int i  = 0  ; i < nstats; i++ ) { 
        double f4ratio = f4ratiodenom[i]>0 ? (f4rationum[i]/f4ratiodenom[i]):0;
        vector<double> &jnum = jf4rationum[i];
        vector<double> &jdenom =   jf4ratiodenom[i];
        vector<double> &jw = jf4ratiow[i];
        vector<double> jf4ratio (jnum.size());
        if (io::debug >= 1)
            cout << "weight = " << w[i] << endl;
        for (int j = 0 ; j < jnum.size(); j++) {
            jnum[j] = f4rationum[i] - jnum[j];
            jdenom[j] = f4ratiodenom[i] - jdenom[j];
            jf4ratio[j] = jdenom[j]>0?(jnum[j]/jdenom[j]):0;
            jw[j] = w[i] - jw[j];
        }
        if (io::debug >= 1) {
            for (int j = 0 ; j < jnum.size(); j++){
                cout << jf4ratio[j] << "\t" << jnum[j] << "\t" << jdenom[j] << "\t" << jw[j] << endl;
            }
        }
        pair<double, double> p = functions::weightedjack ( jf4ratio, jw, f4ratio);
        cout << "result:\t";
        vectorfn::printvector(pops[i]); cout << "\t" ; 
        cout << f4ratio << "\t";
        cout << p.first << "\t" << p.second << "\t" << p.first/p.second << "\t";
//        cout << baba[i] << "\t" << abba[i] << "\t";
        cout << round(baba[i]) << "\t" << round(abba[i]) << "\t" << w[i]  << endl;
//        cout << f4rationum[i] << "\t" << f4ratiodenom[i] << "\t" ; cout << p.first << "\t" << p.second << endl;
    }
}

int main (int argc, char* argv[]) {
	f4ratio f = f4ratio(argc,argv);
    f.run ();
    f.summarize();
}

