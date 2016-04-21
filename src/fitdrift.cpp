#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "vectorfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "convertf.h"
#include "expfit.h"

ofstream statsfs;

class fitdrift  { 
    public:
        fitdrift (int argc, char *argv[]) ;
        void run_theory (double f, double tau) ;
        void read_afs ( string filename , vector<double> &x, vector<double> &y)  ;
        void fit ()  ;

        int col;
        double l ;
        double h;
        double tau;
        double f;
        double taumin  ;
        double taumax ;
        int taubins ; 
        double fmin ;
        double fmax ;
        int fbins;
        int n ;
        int samplesize ; 
        int nmin, nmax;

        vector<double> datax;
        vector<double> datay;
        vector<double> theoryx;
        vector<double> theoryy;
        vector<double> expectedy;


        vector<double> tgrid;
        vector<double> fgrid;

        data *d, *p;
        expfit *e;
        double datalambda;
        double dataintercept;
        double datasigma2;
};

fitdrift::fitdrift (int argc, char *argv[]) {
	static const char *optString = "vh";

	static const struct option longOpts[] = {
		{ "parameter", required_argument, NULL, 'p' },
		{ "verbose", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, no_argument, NULL, 0 }
	};

    d = new data (argc, argv, optString, longOpts);
	string pfile;
	d->get_string ("parameter",pfile,true);
	p = new data (pfile);
	p->print_parameters();
    p->get_int("col",col,1);
    p->get_double("l",l,1);
    p->get_double("h",h,1000);
    p->get_int ("debug", io::debug, 0);
    
    string datafile ;
    p->get_string ("datafile",datafile, true);
    read_afs ( datafile, datax, datay);
    expectedy.resize (datay.size());

    p->get_double ("tau", tau, 0.1);
    p->get_double ("f", f, 0.04);
    e = new expfit (col,l,h,100);
    e->set_data (datax, datay);
    cout << "nmin = " << nmin <<"\tnmax = " << nmax <<endl;
    e->findmle (dataintercept, datalambda, datasigma2);
    cout << "observed mle = " << dataintercept <<","<<datalambda<<endl;

    tgrid = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.20};
//    tgrid = {0.01,0.03,0.05,0.07,0.1,0.2};
    fgrid = {0.01,0.03};

}
 
void fitdrift::run_theory (double f, double tau) { 
    ostringstream oss ;
    oss << "~np29/bin/kimtrans " << f <<" " << tau << ">theory.f-"<<f<<".tau-"<<tau;

    ostringstream ss;
    ss << "theory.f-" <<f<<".tau-"<<tau;

    system (oss.str().c_str());
    ifstream inp (ss.str());
    if (!inp.is_open()){
        cerr << "Error reading file theory"  <<endl;
        exit(1);
    }
    string line;
    theoryx.resize(0);
    theoryy.resize(0);
    double norm ;
    double total   = 0;
    int linenum = 0 ;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#') {
            if (linenum == 2) { 
                vector<string> toks;
                functions::tokenize(line.c_str(),toks," \t");
                double f0 = atof (toks[toks.size()-2].c_str());
                double f1 = atof (toks[toks.size()-1].c_str());
                norm = 1 -f0 -f1;
            }
            continue;
        }
		vector<string> toks;
		functions::tokenize(line.c_str(),toks," \t");
        theoryx.push_back (atof(toks[0].c_str()));
        double y= atof(toks[1].c_str());
        theoryy.push_back(y);
        total += y;
    }
    norm /=total;
    for (int i = 0 ; i < theoryy.size(); i++) 
        theoryy[i] *= norm;

    if (io::debug>=1) {
        for (int i = 0 ; i < theoryy.size(); i++)  {
            cout << theoryx[i] << "\t" << theoryy[i] << endl;
        }
    }

    for (int i = 0 ; i < datax.size();  i++) {  
        double w =0 ;
        int b = datax[i];
        for (int j = 0 ; j < theoryx.size(); j++) {  
            double tmp = theoryx[j];
            w += pow(tmp,b)*pow(1-tmp,samplesize-b)*theoryy[j];
        }
        w *= mathfn::choose (samplesize,b);
        expectedy[i] = w;
    }

    for (int i = 0 ; i < datax.size(); i++) { 
        cout << datax[i] <<"\t"<<datay[i] << "\t"<<expectedy[i] << endl;
    }
    inp.close ();
}

void fitdrift::read_afs ( string filename , vector<double> &x, vector<double> &y)  {
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
    nmin = -1;
    nmax = -1;
    n = 0 ;
    
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		vector<string> toks;
		functions::tokenize(line.c_str(),toks," \t");
        x.push_back (atof(toks[0].c_str()));
        y.push_back (atof(toks[col].c_str()));
        if (x[n]>=l && nmin <0)
            nmin = n;
        if (x[n]>=h && nmax<0)
            nmax = n;
        n++;
	 
    }
    if (nmax  < 0)
        nmax = n - 1;
    samplesize = n  - 1 ;
    cout << "samplesize = " << samplesize << endl;
    inp.close ();
}

void fitdrift::fit ()  {
    vector< vector<double>  > ey; 
    for (int i = 0 ; i < tgrid.size()*fgrid.size(); i++) {
        ey.push_back(vector<double> (datax.size()));
    }

    ofstream ofs ("expected.afs");
    string s1 = "#";
    string s2 = "#";
    for (int i =  0 ;  i < tgrid.size(); i++) {
        for (int k = 0 ;  k < fgrid.size(); k++){
            s1 +=  tostring( tgrid[i]) + "\t" ;
            s2 +=  tostring( fgrid[k]) + "\t" ;
        }
    }
    s1 += "\n";
    s2 += "\n";
    ofs << s1;
    ofs << s2;

    for (int i =  0 ;  i < tgrid.size(); i++) {
        for (int k = 0 ;  k < fgrid.size(); k++){
            run_theory (fgrid[k],tgrid[i]);
            e->set_data (datax, expectedy);
            for (int j = 0 ;  j < expectedy.size(); j++) {
                ey[i][j] = expectedy[j];
            }
            double intercept, lambda, sigma2;
            e->findmle (intercept, lambda, sigma2);
            double diff =  pow( (lambda-datalambda),2);
            cout << "t = " << tgrid[i] << "\tf = " <<fgrid[k] << "\tmle = " << intercept <<","<<lambda<<"\t diff =  "<<diff << endl;

        }
    }
    for (int i = 0 ; i < datax.size(); i++) {
        ofs << datax[i] << "\t";
        for (int j = 0 ; j < tgrid.size()*fgrid.size(); j++)
            ofs << ey[j][i] << "\t";
        ofs << endl;
    }
    ofs.close ();
}




int main (int argc, char *argv[]) { 
    fitdrift fd = fitdrift (argc, argv);
    fd.fit ();
}
