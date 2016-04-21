#include "ram.h"
#include "vectorfn.h"

class computefst { 

    public:
        computefst  (int argc, char *argv[]) ;
    //    ~computefst() ;
        void process ( string pop1, string pop2 );

    vector<int> start; vector<int> end; vector<int> val;
    data *d;
    data *p;
    int window;
    vector<double> fsta;
    vector<double> fstx;
    vector<double> diva;
    vector<double> divx;
    vector<double> cova;
    vector<double> covx;

    string goodsnpname;
    bool givengoodsnps;
    unordered_map <string, string> * gsmap  ;

    bool onlycommon;
    double freqthreshold ; 

    string pop1chra;
    string pop2chra;
    int bin; 

};
/*
computefst::~computefst() { 
    delete p;
    delete d;
}*/

computefst::computefst  (int argc, char *argv[]) { 
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

    p->get_string ("pop1chr", pop1chra, true);
    p->get_string ("pop2chr", pop2chra, true);
    p->get_int ("debug", io::debug, 0);
    p->get_int ("window", window, 10000, false); window = window * 1000;
    p->get_double ("freqthreshold", freqthreshold, -1, false);
	givengoodsnps = p->get_string ("goodsnpname",goodsnpname,false);
    if (freqthreshold > 0)
        onlycommon =  1;
    cout << "debug =  " <<io::debug << endl;

    gsmap  = NULL;
    givengoodsnps  = false;

    if (functions::fileExists (goodsnpname) ) { 
        gsmap =new unordered_map <string, string> ();
        fileio::read_map (goodsnpname, gsmap, 0);
        givengoodsnps  = true;
        cout << "goodsnpname =  " << goodsnpname  << endl;
    }

    process ( pop1chra, pop2chra);

}




void computefst::process ( string pop1, string pop2 ){

    vector<string> p1;
    vector<string> p2;
    vector<double> fstnumv;
    vector<double> fstdenomv;
    vector<double> fstdenomv2;
    vector<double> numsnpsv;

    functions::tokenize (pop1, p1, ",");
    functions::tokenize (pop2, p2, ",");
    if ( !(p1.size()==p2.size()) ) { 
        cerr << "List of files " << pop1 << "," << pop2  <<  " must be equal " << endl;
        exit (1);
    }
    double fstnum = 0 ;
    double fstdenom = 0 ;
    double fstdenom2 = 0 ;
    double tmpfstnum = 0 ;
    double tmpfstdenom = 0 ;
    double tmpfstdenom2 = 0 ;
    int tmpnumsnp =  0;
    double div1 = 0 ; 
    double div2 = 0 ;
    int numsnp = 0 ;

    for (int i  = 0 ;  i < p1.size(); i ++)  {
        if (io::debug>=0)
            cout << "Opening file   " << p1[i] << endl;
        ifstream p1afs (p1[i].c_str());
        ifstream p2afs (p2[i].c_str());

        string line1;
        string line2;
        bool y1 = std::getline (p1afs,line1).good();
        bool y2 = std::getline (p2afs,line2).good();
        int wnum  = -1;

       int linenum = 0 ;
 
        while (y1 && y2 ) { 
            linenum ++;
            if (io::debug>=1)
                cout << "Linenum = " << linenum << endl;
            bool flag = true;
            while (y1 && line1[0]=='#') { 
                y1 = std::getline (p1afs,line1).good();
            }
            while (y2 && line2[0]=='#') { 
                y2 = std::getline (p2afs,line1).good();
            }
            vector<string> toks1;
            functions::tokenize (line1.c_str(), toks1, " \t");
            string id1 = toks1[0];
            int a1 = atoi(toks1[1].c_str()); 
            int b1 = atoi(toks1[2].c_str());
            int n1 = a1 + b1;
            double f1 = (1.*a1)/n1;

            vector<string> toks2;
            functions::tokenize (id1.c_str(), toks2, ":");
            string chr = toks2[0];
            long pos = atol (toks2[1].c_str());

            toks1.resize(0);
            functions::tokenize (line2.c_str(), toks1, " \t");
            string id2  = toks1[0].c_str();
            int a2 = atoi(toks1[1].c_str()); 
            int b2 = atoi(toks1[2].c_str());
            int n2 = a2 + b2;
            double f2 = (1.*a2)/n2;

           

            if ( !(id1.compare(id2)==0 )) {
                cerr << "Files " << p1[i] << "," << p2[i] <<  " differ at line " << linenum << endl;
                exit (1);
            } else  {
            }

            if (givengoodsnps ){
                if (gsmap->find(id1)==gsmap->end())
                    flag = false;
            }

            if (io::debug >= 1) 
                cout << "There " << flag  << "\t" << wnum  << endl;

            if (wnum == -1 || (wnum + 1)  * window <= pos ) { 
                if (wnum >=0 ) {

                    if (tmpfstdenom > 0 ) { 
                        double localfst = tmpfstdenom > 0? tmpfstnum/tmpfstdenom:0;
                        fstnumv.push_back (tmpfstnum);
                        fstdenomv.push_back (tmpfstdenom);
                        numsnpsv.push_back(tmpnumsnp);
                    }
                }
                wnum = pos/window;
                tmpfstnum = 0 ; tmpfstdenom = 0 ;
                tmpnumsnp = 0 ;
            }
 
            if (flag ) { 
                double h1 = (1.*a1*(n1-a1))/(n1*(n1-1)); 
                double h2 = (1.*a2*(n2-a2))/(n2*(n2-1)); 
                double tmp = pow( ((1.*a1)/n1-(1.*a2)/n2),2) - h1/n1 - h2/n2;
                if (!onlycommon ||( f1>=freqthreshold && (1-f1) >= freqthreshold && f2 >=freqthreshold && (1-f2) >= freqthreshold)){
//                if ( f1 >=0.05 && f1<= 0.95 && f2 >=0.05 && f2 <=0.95){
                    fstnum += tmp;
                    fstdenom += (tmp + h1 + h2);
                    fstdenom2 += 0.5*tmp + h1 + h2 ;

                    tmpfstnum += tmp;
                    tmpfstdenom += (tmp + h1 + h2);
                }
                numsnp ++;
                tmpnumsnp ++;
            }

            y1 = std::getline (p1afs,line1).good();
            y2 = std::getline (p2afs,line2).good();
        }
        p1afs.close ();
        p2afs.close ();
    }
    double fst = fstnum/fstdenom;
    double fst2 = fstnum/fstdenom2;
    vector<double> fstv;
    for (int i =  0 ; i < fstnumv.size(); i++) { 
        fstnumv[i] = fstnum - fstnumv[i];
        fstdenomv[i] = fstdenom - fstdenomv[i];
        numsnpsv[i] = numsnp - numsnpsv[i];
        fstv.push_back(fstnumv[i]/fstdenomv[i]);
    }
    pair<double,double> p = functions::weightedjack (fstv, numsnpsv, fst);

    cout << fst << "\t" << p.first << "\t" << p.second << endl;
    cout << fst2 << endl;
}


int main (int argc, char *argv[]) { 
    computefst c (argc, argv);
}
