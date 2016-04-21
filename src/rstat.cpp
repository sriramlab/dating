#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "printable.h"
#include "genotype.h"
#include "genomecompute.h"


class rstat  : public genomecompute { 
	public:
		rstat  (int argc, char *argv[]) ;
        void compute (vector<int> &gtype, int jackindex) ;
        void setparams () ;
        void summarize ();

	    int nstats;
        vector<double> baba;
        vector<double> abba;
        vector<double> rstatnum ; 
        vector<double> rstatdenom ; 
        vector<double> dw;
        vector< vector<string> > pops;
        vector< vector< double> > jrstatnum;
        vector< vector< double> > jrstatdenom;
        vector< vector< double> > jrstatw;

};



rstat::rstat (int argc, char *argv[]) 
    : genomecompute (argc, argv)  {

}


void rstat::setparams () { 
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

    rstatnum.resize (nstats);
    rstatdenom.resize (nstats);
    dw.resize (nstats);
    jrstatnum.resize ( nstats, vector<double> (njbins) ) ;
    jrstatdenom.resize ( nstats, vector<double> (njbins) ) ;
    jrstatw.resize (nstats, vector<double> (njbins) );
    for (int i  = 0  ; i < nstats; i ++ ) {
        vectorfn::printvector(pops[i]); cout << endl;
    }

}

void rstat::compute (vector<int> &gtype, int jackindex) {
    for (int i  = 0  ; i < nstats; i++ ) { 
        vector<string> &p = pops[i];
        double a = g->get_freq (gtype,  p[0]); 
        double b = g->get_freq (gtype,  p[1]); 
        double c = g->get_freq (gtype,  p[2]); 
        double d = g->get_freq (gtype,  p[3]); 
        double e = g->get_freq (gtype,  p[4]); 
        if (a<0 || b<0 || c<0 || d<0 || e <0 ) {
           continue;        
        }

        double w = (1-c)*d;
        double dnum = (a-b) * (c-d);
        double ddenom = (e-b) * (c-d);
        rstatnum[i] += dnum;
        rstatdenom[i] += ddenom;
        dw[i] += 1;
        jrstatnum[i][jackindex] += dnum;
        jrstatdenom[i][jackindex] += ddenom;
        jrstatw[i][jackindex] += 1;
    }

}


void rstat::summarize () {
    cout << "#Populations";
    cout << "\tR-statistic\tBias-corrected R-statisic\tSE\tZ-score";
    cout << endl;
    for (int i  = 0  ; i < nstats; i++ ) { 
        double rstat = abs(rstatdenom[i])>0 ? (rstatnum[i]/rstatdenom[i]):0;
        vector<double> &jnum = jrstatnum[i];
        vector<double> &jdenom =   jrstatdenom[i];
        vector<double> &jw = jrstatw[i];
        vector<double> jrstat (jnum.size());
        if (io::debug >= 1)
            cout << "weight = " << dw[i] << endl;
        for (int j = 0 ; j < jnum.size(); j++) {
            jnum[j] = rstatnum[i] - jnum[j];
            jdenom[j] = rstatdenom[i] - jdenom[j];
            jrstat[j] = abs(jdenom[j])>0?(jnum[j]/jdenom[j]):0;
            jw[j] = dw[i] - jw[j];
        }
        if (io::debug >= 1) {
            for (int j = 0 ; j < jnum.size(); j++){
                cout << jrstat[j] << "\t" << jnum[j] << "\t" << jdenom[j] << "\t" << jw[j] << endl;
            }
        }
        pair<double, double> p = functions::weightedjack ( jrstat, jw, rstat);
        cout << "result:\t";
        vectorfn::printvector(pops[i]); cout << "\t" ; 
        cout << rstat << "\t";
        cout << p.first << "\t" << p.second << "\t" << p.first/p.second << "\t";
//        cout << baba[i] << "\t" << abba[i] << "\t";
//        cout << round(baba[i]) << "\t" << round(abba[i]) << "\t" << dw[i]  << endl;
        cout << rstatnum[i] << "\t" << rstatdenom[i] << "\t" ; cout << p.first << "\t" << p.second << endl;
    }
}

int main (int argc, char* argv[]) {
	rstat d = rstat(argc,argv);
    d.run ();
    d. summarize();
}

