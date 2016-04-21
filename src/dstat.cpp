#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "dstat.h"

dstat::dstat (int argc, char *argv[]) 
    : genomecompute (argc, argv)  {

}


void dstat::setparams () { 
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

    dstatnum.resize (nstats);
    dstatdenom.resize (nstats);
    dw.resize (nstats);
    jdstatnum.resize ( nstats, vector<double> (njbins) ) ;
    jdstatdenom.resize ( nstats, vector<double> (njbins) ) ;
    jdstatw.resize (nstats, vector<double> (njbins) );
    for (int i  = 0  ; i < nstats; i ++ ) {
        vectorfn::printvector(pops[i]); cout << endl;
    }

}

void dstat::compute (vector<int> &gtype, int jackindex) {
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
        dstatnum[i] += dnum;
        dstatdenom[i] += ddenom;
        dw[i] += 1;
        jdstatnum[i][jackindex] += dnum;
        jdstatdenom[i][jackindex] += ddenom;
        jdstatw[i][jackindex] += 1;
        if (io::debug>=1) { 
            cout << baba[i] << "\t" << abba[i] << endl;
        }
    }

}


void dstat::summarize () {
    cout << "#Populations";
    cout << "\tD-statistic\tBias-corrected D-statisic\tSE\tZ-score";
    cout << "\tABBA\tBABA\tWeights";
    cout << endl;
    for (int i  = 0  ; i < nstats; i++ ) { 
        double dstat = dstatdenom[i]>0 ? (dstatnum[i]/dstatdenom[i]):0;
        vector<double> &jnum = jdstatnum[i];
        vector<double> &jdenom =   jdstatdenom[i];
        vector<double> &jw = jdstatw[i];
        vector<double> jdstat (jnum.size());
        if (io::debug >= 1)
            cout << "weight = " << dw[i] << endl;
        for (int j = 0 ; j < jnum.size(); j++) {
            jnum[j] = dstatnum[i] - jnum[j];
            jdenom[j] = dstatdenom[i] - jdenom[j];
            jdstat[j] = jdenom[j]>0?(jnum[j]/jdenom[j]):0;
            jw[j] = dw[i] - jw[j];
        }
        if (io::debug >= 1) {
            for (int j = 0 ; j < jnum.size(); j++){
                cout << jdstat[j] << "\t" << jnum[j] << "\t" << jdenom[j] << "\t" << jw[j] << endl;
            }
        }
        pair<double, double> p = functions::weightedjack ( jdstat, jw, dstat);
        cout << "result:\t";
        vectorfn::printvector(pops[i]); cout << "\t" ; 
        cout << dstat << "\t";
        cout << p.first << "\t" << p.second << "\t" << p.first/p.second << "\t";
//        cout << baba[i] << "\t" << abba[i] << "\t";
        cout << round(baba[i]) << "\t" << round(abba[i]) << "\t" << dw[i]  << endl;
//        cout << dstatnum[i] << "\t" << dstatdenom[i] << "\t" ; cout << p.first << "\t" << p.second << endl;
    }
}

int main (int argc, char* argv[]) {
	dstat d = dstat(argc,argv);
    d.run ();
    d.summarize();
}

