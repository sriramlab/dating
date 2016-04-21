#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "indanc.h"

ofstream freqfs;
indanc::indanc (int argc, char *argv[]) 
    : genomecompute (argc, argv)  {

}


void indanc::setparams () { 
    string poplist = "";
    p->get_string ("popfilename", poplist, true);

    freqout = false;
    freqout  = p->get_string ("freqout", freqoutfile, false);
    if (freqout )  {
        freqfs.open (freqoutfile.c_str());
    }

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
        pops.push_back (line);
    }
    inp.close ();
    nstats = pops.size();

    indancnum.resize (nstats);
    indancdenom.resize (nstats);
    w.resize (nstats);
    jindancnum.resize ( nstats, vector<double> (njbins) ) ;
    jindancdenom.resize ( nstats, vector<double> (njbins) ) ;
    jindancw.resize (nstats, vector<double> (njbins) );
    for (int i  = 0  ; i < nstats; i ++ ) {
        cout << pops[i]; cout << endl;
    }

}

void indanc::compute (vector<int> &gtype, int jackindex) {
    for (int i  = 0  ; i < nstats; i++ ) { 
        string p = pops[i];
        double a = g->get_freq (gtype,  p); 


        if (freqout ) { 
            freqfs << a << "\t";
        }

        if (a<0 ) {
           continue;        
        }
        indancnum[i] += a;
        indancdenom[i] ++;
        w[i] ++;

        jindancnum[i][jackindex] += a;
        jindancdenom[i][jackindex] ++;
        jindancw[i][jackindex]++;

    }
    if (freqout ) { 
        freqfs << endl;
    }

}


void indanc::summarize () {
    cout << "#Population";
    cout << "\tAncestry\tBias-corrected ancestry\tSE\tZ-score\tNumerator\tDenominator";
    cout << endl;
    for (int i  = 0  ; i < nstats; i++ ) { 
        double indanc = indancdenom[i]>0 ? (indancnum[i]/indancdenom[i]):0;
        vector<double> &jnum = jindancnum[i];
        vector<double> &jdenom =   jindancdenom[i];
        vector<double> &jw = jindancw[i];
        vector<double> jindanc (jnum.size());
        if (io::debug >= 1)
            cout << "weight = " << w[i] << endl;
        
        for (int j = 0 ; j < jnum.size(); j++) {
            jnum[j] = indancnum[i] - jnum[j];
            jdenom[j] = indancdenom[i] - jdenom[j];
            jindanc[j] = jdenom[j]>0?(jnum[j]/jdenom[j]):0;
        }

        if (outputjack) { 
            string filename = outputjackstring + "." + tostring (i);
            ofstream ofs (filename.c_str());
            for (int j = 0 ; j < jnum.size(); j++) {
                ofs << j << "\t" << jw[j] << "\t" << jindanc[j]<<endl;
            }
            ofs.close ();
        }
        if (io::debug >= 1) {
            for (int j = 0 ; j < jnum.size(); j++){
                cout << jindanc[j] << "\t" << jnum[j] << "\t" << jdenom[j] << "\t" << jw[j] << endl;
            }
        }
        pair<double, double> p = functions::weightedjack ( jindanc, jw, indanc);
        cout << "result:\t";
        cout << pops[i]; cout << "\t" ; 
        cout << indanc << "\t";
        cout << p.first << "\t" << p.second << "\t" << p.first/p.second << "\t" << indancnum[i] << "\t" << indancdenom[i] << "\n";
    }
    if (freqout ) { 
        freqfs.close();
    }
}

int main (int argc, char* argv[]) {
	indanc d = indanc(argc,argv);
    d.run ();
    d.summarize();
}

