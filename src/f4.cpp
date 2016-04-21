#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "f4.h"

f4::f4 (int argc, char *argv[]) 
    : genomecompute (argc, argv)  {

}


void f4::setparams () { 
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

    normpop =  "";
    normalize = p->get_string ("normpop", normpop, false);

    nstats = pops.size();
    abba.resize (nstats);
    baba.resize (nstats);

    f4num.resize (nstats);
    f4denom.resize (nstats);
    dw.resize (nstats);
    jf4num.resize ( nstats, vector<double> (njbins) ) ;
    jf4denom.resize ( nstats, vector<double> (njbins) ) ;
    jf4w.resize (nstats, vector<double> (njbins) );
    for (int i  = 0  ; i < nstats; i ++ ) {
        vectorfn::printvector(pops[i]); cout << endl;
    }

}

void f4::compute (vector<int> &gtype, int jackindex) {
    for (int i  = 0  ; i < nstats; i++ ) { 
        vector<string> &p = pops[i];
        double a = g->get_freq (gtype,  p[0]); 
        double b = g->get_freq (gtype,  p[1]); 
        double c = g->get_freq (gtype,  p[2]); 
        double d = g->get_freq (gtype,  p[3]); 

        if (a<0 || b<0 || c<0 || d<0) {
           continue;        
        }

        double e = 1;
        if ( normalize ) { 
            e = g->get_freq (gtype, normpop); 
            if (e < 0)
               continue;
           e = e * (1-e); 
        }



        double w = (1-c)*d;
        double dnum = (a-b) * (c-d);
        baba[i] += (a*(1-b)*c*(1-d) + (1-a)*b*(1-c)*d);
        abba[i] += ((1-a)*b*c*(1-d) + a*(1-b)*(1-c)*d);
        f4num[i] += dnum;
        f4denom[i] += e;
        dw[i] += 1;
        jf4num[i][jackindex] += dnum;
        jf4denom[i][jackindex] += e;
        jf4w[i][jackindex] += 1;
        if (io::debug>=1) { 
            cout << baba[i] << "\t" << abba[i] << endl;
        }

    }

}


void f4::summarize () {
    cout << "#Populations";
    cout << "\tF4-statistic\tBias-corrected F4-statistic\tSE\tZ-score";
    cout << "\tABBA\tBABA\tWeights";
    cout << endl;
    for (int i  = 0  ; i < nstats; i++ ) { 
        double f4 = f4denom[i]>0 ? (f4num[i]/f4denom[i]):0;
        vector<double> &jnum = jf4num[i];
        vector<double> &jdenom =   jf4denom[i];
        vector<double> &jw = jf4w[i];
        vector<double> jf4 (jnum.size());
        if (io::debug >= 1)
            cout << "weight = " << dw[i] << endl;
        for (int j = 0 ; j < jnum.size(); j++) {
            jnum[j] = f4num[i] - jnum[j];
            jdenom[j] = f4denom[i] - jdenom[j];
            jf4[j] = jdenom[j]>0?(jnum[j]/jdenom[j]):0;
            jw[j] = dw[i] - jw[j];
        }
        if (outputjack) { 
            string filename = outputjackstring + "." + tostring (i);
            ofstream ofs (filename.c_str());
            for (int j = 0 ; j < jnum.size(); j++) {
                ofs << j << "\t" << jw[j] << "\t" << jf4[j]<<endl;
            }
            ofs.close ();
        }
        if (io::debug >= 1) {
            for (int j = 0 ; j < jnum.size(); j++){
                cout << jf4[j] << "\t" << jnum[j] << "\t" << jdenom[j] << "\t" << jw[j] << endl;
            }
        }
        pair<double, double> p = functions::weightedjack ( jf4, jw, f4);
        cout << "result:\t";
        vectorfn::printvector(pops[i]); cout << "\t" ; 
        cout << f4 << "\t";
        cout << p.first << "\t" << p.second << "\t" << p.first/p.second << "\t";
//        cout << baba[i] << "\t" << abba[i] << "\t";
        cout << round(baba[i]) << "\t" << round(abba[i]) << "\t" << dw[i]  << endl;
//        cout << f4num[i] << "\t" << f4denom[i] << "\t" ; cout << p.first << "\t" << p.second << endl;
    }
}

int main (int argc, char* argv[]) {
	f4 d = f4(argc,argv);
    d.run ();
    d.summarize();
}

