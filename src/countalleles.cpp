#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "countalleles.h"

countalleles::countalleles (int argc, char *argv[]) 
    : genomecompute (argc, argv)  {

}


void countalleles::setparams () { 
    countallelesnum.resize (nstats);
    countallelesdenom.resize (nstats);
    w.resize (nstats);
    jcountallelesnum.resize ( nstats, vector<double> (njbins) ) ;
    jcountallelesdenom.resize ( nstats, vector<double> (njbins) ) ;
    jcountallelesw.resize (nstats, vector<double> (njbins) );
    for (int i  = 0  ; i < nstats; i ++ ) {
        cout << pops[i]; cout << endl;
    }

}

void countalleles::compute (vector<int> &gtype, int jackindex) {
    for (int i  = 0  ; i < nstats; i++ ) { 
        string p = pops[i];
        double a = g->get_freq (gtype,  p); 
        if (a<0 ) {
           continue;        
        }
        countallelesnum[i] += a;
        countallelesdenom[i] ++;
        w[i] ++;

        jcountallelesnum[i][jackindex] += a;
        jcountallelesdenom[i][jackindex] ++;
        jcountallelesw[i][jackindex]++;

    }

}


void countalleles::summarize () {
    cout << "#Population";
    cout << "\tAncestry\tBias-corrected ancestry\tSE\tZ-score";
    cout << endl;
    for (int i  = 0  ; i < nstats; i++ ) { 
        double countalleles = countallelesdenom[i]>0 ? (countallelesnum[i]/countallelesdenom[i]):0;
        vector<double> &jnum = jcountallelesnum[i];
        vector<double> &jdenom =   jcountallelesdenom[i];
        vector<double> &jw = jcountallelesw[i];
        vector<double> jcountalleles (jnum.size());
        if (io::debug >= 1)
            cout << "weight = " << w[i] << endl;
        
        for (int j = 0 ; j < jnum.size(); j++) {
            jnum[j] = countallelesnum[i] - jnum[j];
            jdenom[j] = countallelesdenom[i] - jdenom[j];
            jcountalleles[j] = jdenom[j]>0?(jnum[j]/jdenom[j]):0;
        }

        if (outputjack) { 
            string filename = outputjackstring + "." + tostring (i);
            ofstream ofs (filename.c_str());
            for (int j = 0 ; j < jnum.size(); j++) {
                ofs << j << "\t" << jw[j] << "\t" << jcountalleles[j]<<endl;
            }
            ofs.close ();
        }
        if (io::debug >= 1) {
            for (int j = 0 ; j < jnum.size(); j++){
                cout << jcountalleles[j] << "\t" << jnum[j] << "\t" << jdenom[j] << "\t" << jw[j] << endl;
            }
        }
        pair<double, double> p = functions::weightedjack ( jcountalleles, jw, countalleles);
        cout << "result:\t";
        cout << pops[i]; cout << "\t" ; 
        cout << countalleles << "\t";
        cout << p.first << "\t" << p.second << "\t" << p.first/p.second << "\n";
    }
}

int main (int argc, char* argv[]) {
	countalleles d = countalleles(argc,argv);
    d.run ();
    d.summarize();
}

