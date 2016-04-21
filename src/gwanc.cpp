#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "gwanc.h"

gwanc::gwanc (int argc, char *argv[]) 
    : genomecompute (argc, argv)  {

}


void gwanc::setparams () { 
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

    gwancnum.resize (nstats);
    gwancdenom.resize (nstats);
    dw.resize (nstats);
    jgwancnum.resize ( nstats, vector<double> (njbins) ) ;
    jgwancdenom.resize ( nstats, vector<double> (njbins) ) ;
    jgwancw.resize (nstats, vector<double> (njbins) );
    for (int i  = 0  ; i < nstats; i ++ ) {
        vectorfn::printvector(pops[i]); cout << endl;
    }

}

void gwanc::compute (vector<int> &gtype, int jackindex) {
    for (int i  = 0  ; i < nstats; i++ ) { 
        vector<string> &p = pops[i];
        double a = g->get_freq (gtype,  p[0]); 
        double b = g->get_freq (gtype,  p[1]); 
        double c = g->get_freq (gtype,  p[2]); 
        double d = g->get_freq (gtype,  p[3]); 
        if (a<0 || b<0 || c<0 || d<0) {
           continue;        
        }

        double w = (1-c)*d;
        double dnum = (a-b) * (c-d);
        double ddenom = (a+b-2*a*b) * (c+d-2*c*d);
        baba[i] += (a*(1-b)*c*(1-d) + (1-a)*b*(1-c)*d);
        abba[i] += ((1-a)*b*c*(1-d) + a*(1-b)*(1-c)*d);
        gwancnum[i] += dnum;
        gwancdenom[i] += ddenom;
        dw[i] += 1;
        jgwancnum[i][jackindex] += dnum;
        jgwancdenom[i][jackindex] += ddenom;
        jgwancw[i][jackindex] += 1;
    }

}


void gwanc::summarize () {
    cout << "#Populations";
    cout << "\tD-statistic\tBias-corrected D-statisic\tSE\tZ-score";
    cout << "\tABBA\tBABA\tWeights";
    cout << endl;
    for (int i  = 0  ; i < nstats; i++ ) { 
        double gwanc = gwancdenom[i]>0 ? (gwancnum[i]/gwancdenom[i]):0;
        vector<double> &jnum = jgwancnum[i];
        vector<double> &jdenom =   jgwancdenom[i];
        vector<double> &jw = jgwancw[i];
        vector<double> jgwanc (jnum.size());
        if (io::debug >= 1)
            cout << "weight = " << dw[i] << endl;
        for (int j = 0 ; j < jnum.size(); j++) {
            jnum[j] = gwancnum[i] - jnum[j];
            jdenom[j] = gwancdenom[i] - jdenom[j];
            jgwanc[j] = jdenom[j]>0?(jnum[j]/jdenom[j]):0;
            jw[j] = dw[i] - jw[j];
        }
        if (io::debug >= 1) {
            for (int j = 0 ; j < jnum.size(); j++){
                cout << jgwanc[j] << "\t" << jnum[j] << "\t" << jdenom[j] << "\t" << jw[j] << endl;
            }
        }
        pair<double, double> p = functions::weightedjack ( jgwanc, jw, gwanc);
        cout << "result:\t";
        vectorfn::printvector(pops[i]); cout << "\t" ; 
        cout << gwanc << "\t";
        cout << p.first << "\t" << p.second << "\t" << p.first/p.second << "\t";
//        cout << baba[i] << "\t" << abba[i] << "\t";
        cout << round(baba[i]) << "\t" << round(abba[i]) << "\t" << dw[i]  << endl;
//        cout << gwancnum[i] << "\t" << gwancdenom[i] << "\t" ; cout << p.first << "\t" << p.second << endl;
    }
}

int main (int argc, char* argv[]) {
	gwanc d = gwanc(argc,argv);
    d.run ();
    d.summarize();
}

