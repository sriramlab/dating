#include "ram.h"
#include "vectorfn.h"

class computexaratio { 

    public:
        computexaratio  (int argc, char *argv[]) ;
    //    ~computexaratio() ;
        void process ( string pop1, string pop2, string covariates, string output, vector<double> &fst, vector<double> &div, vector<double> &cov  );
        void process ( string pop1, string pop2, string pop3, string covariates, string output, vector<double> &drift12, vector<double> &drift23, vector<double> &drift31, vector<double> &div, vector<double> &cov  );
        void summary(int);

    vector<int> start; vector<int> end; vector<int> val;
    data *d;
    data *p;
    int window;
    double gwindow;
    vector<double> fsta;
    vector<double> fstx;
    vector<double> diva;
    vector<double> divx;
    vector<double> cova;
    vector<double> covx;

    vector<double> adrift1, adrift2, adrift3;
    vector<double> xdrift1, xdrift2, xdrift3;
    bool givenpop3;

    string pop1chra;
    string pop2chra;
    string pop3chra;
    string pop1chrx;
    string pop2chrx;
    string pop3chrx;
    string covariatechra;
    string covariatechrx;
    string outputa;
    string outputx;
    string prefix;
    int bin; 
    double covamin, covamax;
    double covxmin, covxmax;
    double covmin, covmax;


    bool givensnps;
    bool usegenpos;
    string snpfile;
    unordered_map <string, string> *snpmap;

};
/*
computexaratio::~computexaratio() { 
    delete p;
    delete d;
}*/

computexaratio::computexaratio  (int argc, char *argv[]) { 
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

    p->get_string ("pop1chra", pop1chra, true);
    p->get_string ("pop2chra", pop2chra, true);
    givenpop3 = p->get_string ("pop3chra", pop3chra, false);

    p->get_string ("pop1chrx", pop1chrx, true);
    p->get_string ("pop2chrx", pop2chrx, true);
    p->get_string ("pop3chrx", pop3chrx, givenpop3);
    
    p->get_string ("covariatechra", covariatechra, true);
    p->get_string ("covariatechrx", covariatechrx, true);
    p->get_int ("window", window, true);
    p->get_double ("gwindow", gwindow, false);
    p->get_string ("outputa", outputa, true);
    p->get_string ("outputx", outputx, true);
    p->get_string ("prefix", prefix, true);

    givensnps = false;
    givensnps = p->get_string ("snpfile", snpfile, false);
    if (givensnps) { 
        snpmap  = new unordered_map <string, string> ();
        fileio::read_map (snpfile, snpmap, 0);
    }
    p->get_boolean ("usegenpos", usegenpos, false,false);

    p->get_int ("bin", bin, 10, false);
    p->get_int ("debug", io::debug, 0);

    if (usegenpos)
        gwindow = gwindow/100;
    else
        window = window * 1000;

    cout << "debug =  " <<io::debug << endl;

    if (givenpop3) {
        process ( pop1chra, pop2chra, pop3chra, covariatechra, outputa, adrift1, adrift2, adrift3, diva, cova);
        process ( pop1chrx, pop2chrx, pop3chrx, covariatechrx, outputx, xdrift1, xdrift2, xdrift3, divx, covx);
    } else {
        process ( pop1chra, pop2chra, covariatechra, outputa, fsta, diva, cova);
        process ( pop1chrx, pop2chrx, covariatechrx, outputx, fstx, divx, covx);

        pair<double,double> tmp;
        tmp = vectorfn::min(cova); covamin = tmp.first;
        tmp = vectorfn::max(cova); covamax = tmp.first;
        tmp = vectorfn::min(covx); covxmin = tmp.first;
        tmp = vectorfn::max(covx); covxmax = tmp.first;
        covmin = covamin<covxmin ? covamin:covxmin; 
        covmax = covamax>covxmax ? covamax:covxmax; 

        summary(0);
        summary(1);
        summary(2);
        summary(3);
    }
}

void computexaratio::summary (int mode)  {
    vector<double> meancova (bin);
    vector<double> meancovx (bin);
    vector<double> meanfsta (bin);
    vector<double> meanfstx (bin);
    vector<double> meandiva (bin);
    vector<double> meandivx (bin);
    vector<int> denoma (bin);
    vector<int> denomx (bin);


    vector<double> lima (bin-1);
    vector<double> limx (bin-1);
    switch (mode) { 
        case 1:
            {
                double delta  = (covmax-covmin)/bin;
                for (int j  = 0 ; j < bin-1; j++)  {
                    lima[j] = (j+1)*delta; limx[j] = lima[j];
                }
            }
        break;
        case 2:
        {
            vector<double> q (bin-1);

            for (int  j =  0 ; j < bin - 1; j++) 
                q[j] = (1.*(j+1))/bin;
            vectorfn::quantile (covx, q, limx);
            for (int j  = 0 ; j < bin-1; j++)  
                lima[j] = limx[j];
        }
        break;
        case 3:
        {
            vector<double> q (bin-1);

            for (int  j =  0 ; j < bin - 1; j++) 
                q[j] = (1.*(j+1))/bin;
            vectorfn::quantile (covx, q, limx);
            vectorfn::quantile (cova, q, lima);
        }
        break;

        default:
        {
            for (int j  = 0 ; j < bin-1; j++)  {
                limx[j] = (j+1)*100; lima[j] = limx[j];
            }
            break;
        }

    }

    for (int i = 0 ;  i < cova.size() ; i++) { 
        int index =  0;
        for (;index < lima.size() && cova[i] > lima[index] ; index++);
//        cout << "index = " << index <<endl;

        meancova[index] += cova[i];
        meanfsta[index] += fsta[i];
        meandiva[index] += diva[i];
        denoma [index]++;
    }

    for (int i = 0 ;  i < covx.size() ; i++) { 
        int index =  0;
        for (;index < limx.size() && covx[i] > limx[index] ; index++);
//        cout << "index = " << index <<endl;

        meancovx[index] += covx[i];
        meanfstx[index] += fstx[i];
        meandivx[index] += divx[i];
        denomx [index] ++;
    }



    cout << "***Mode " << mode << "***" << endl;
    ostringstream oss ; oss << prefix << ".mode-" << mode << ".bins";

    ofstream ofs ( oss.str());
    ofs.setf(ios::fixed,ios::floatfield); 
	ofs.precision(3);
    for (int i = 0 ;   i < bin ; i++){ 
        double la = 0 ; 
        double lx = 0 ;
        double ra = covamax ;
        double rx = covxmax ;
        if ( i == 0) {
            ra = lima[i]; rx = limx[i];
        } else if (i==bin-1) { 
            la = lima[i-1]; lx = limx[i-1];
        } else {
            la = lima[i-1]; lx = limx[i-1];
            ra = lima[i]; rx = limx[i];
        }

        if (denoma[i] > 0) {
            meancova[i] /= denoma[i]; meanfsta[i] /= denoma[i]; meandiva[i] /= denoma[i];
        }
        if (denomx[i] > 0) {
            meancovx[i] /= denomx[i]; meanfstx[i] /= denomx[i]; meandivx[i] /= denomx[i];
        }
        ofs << i << "\t" << la << "\t" << ra << "\t" << lx << "\t" << rx ;
        ofs << "\t" << meancova[i] << "\t" << meancovx[i] << "\t" << meanfsta[i] << "\t" << meanfstx[i] << "\t" << log(1-meanfsta[i])/log(1-meanfstx[i]) ;
        ofs << "\t" << meandiva[i] << "\t" << meandivx[i] << "\t" << meandivx[i]/meandiva[i]  ;
        ofs << endl;
    }
    ofs.close ();
}


void computexaratio::process ( string pop1, string pop2, string covariates, string output, vector<double> &fst, vector<double> &div, vector<double> &cov ){

    vector<string> p1;
    vector<string> p2;
    vector<string> c;

    functions::tokenize (pop1, p1, ",");
    functions::tokenize (pop2, p2, ",");
    functions::tokenize (covariates, c, ",");
    if ( !(p1.size()==p2.size() && p1.size()==c.size())) { 
        cerr << "List of files " << pop1 << "," << pop2  << "," << covariates << " must be equal " << endl;
        exit (1);
    }
    ofstream ofs (output.c_str());
    ofs.setf(ios::fixed,ios::floatfield); 
	ofs.precision(3);

    for (int i  = 0 ;  i < p1.size(); i ++)  {
        if (io::debug>=0)
            cout << "Opening file   " << p1[i] << endl;
        ifstream p1afs (p1[i].c_str());
        ifstream p2afs (p2[i].c_str());
        ifstream cafs (c[i].c_str());


        string line1;
        string line2;
        string line3;
        string line4;

        bool y1 = std::getline (p1afs,line1).good();
        bool y2 = std::getline (p2afs,line2).good();
        bool y3 = std::getline (cafs,line3).good();


        int wnum  = -1;

        double fstnum = 0 ;
        double fstdenom = 0 ;


        double div1 = 0 ; 
        double div2 = 0 ;
        double avgcov =  0;
        int numsnp = 0 ;
        int linenum = 0 ;
 
        while (y1 && y2 && y3) { 

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
            while (y3 && line3[0]=='#') { 
                y3 = std::getline (cafs,line1).good();
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

            double genpos;
            if (usegenpos) { 
                if (snpmap->find(id1)!=snpmap->end ()){
                    string l = (*snpmap)[id1];
                    vector<string> ltoks;
                    functions::tokenize (l.c_str(), ltoks, " \t");
                    genpos = atof(ltoks[2].c_str());
                } else {
                    cerr << "Cannot find SNP " << id1 <<endl;
                }        
            }

            toks1.resize(0);
            functions::tokenize (line2.c_str(), toks1, " \t");
            string id2  = toks1[0].c_str();
            int a2 = atoi(toks1[1].c_str()); 
            int b2 = atoi(toks1[2].c_str());
            int n2 = a2 + b2;
            double f2 = (1.*a2)/n2;

            toks1.resize(0);
            functions::tokenize (line3.c_str(), toks1, " \t");
            string id3  = toks1[0];
            double val =  atof(toks1[1].c_str());
            if (io::debug >= 1) 
                cout << "Here " << id1 << "\t" << id3 << endl;
            


            if ( !(id1.compare(id2)==0 )) {
                cerr << "Files " << p1[i] << "," << p2[i] << "," << c[i] << " differ at line " << linenum << endl;
                exit (1);
            } else  {
                if ( !(id1.compare (id3)==0)){
                    toks2.resize(0);
                    functions::tokenize (id3.c_str(), toks2, ":");
                    long pos3 = atol (toks2[1].c_str());
                    if ( pos3 < pos )  { 
                        y3 = std::getline (cafs,line3).good();
                        continue;
                    } else 
                        flag = false;
                }
            }

            if (io::debug >= 1) 
                cout << "There " << flag  << "\t" << wnum  << "\t" << genpos <<  endl;


            if (usegenpos ) {
                if (wnum == -1 || (wnum + 1)  * gwindow <= genpos ) { 
                    if (wnum >=0 && numsnp > 10) {
                        double localfst = fstdenom > 0? fstnum/fstdenom:0;
                        double divratio = div2>0 ?div1/div2:0;
                        avgcov /= numsnp;
                        ofs << id1 << "\t"<< chr << "\t" << ((double)(100*gwindow*wnum)) << "\t" << ((double)(100*gwindow*(wnum+1))) << "\t" << ((double)(100*gwindow*(wnum+0.5))) << "\t" << numsnp << "\t" << localfst << "\t" << "\t" << div1 << "\t" << div2 << "\t" <<divratio << "\t" << avgcov << endl;
                        fst.push_back (localfst);
                        div.push_back(divratio);
                        cov.push_back (avgcov);
                    }
                    wnum = genpos/gwindow;
                    fstnum = 0 ; fstdenom = 0 ;
                    div1 = 0 ; div2 = 0 ;
                    avgcov =  0;
                    numsnp = 0 ;
                }

            } else {
                if (wnum == -1 || (wnum + 1)  * window <= pos ) { 
                    if (wnum >=0 && numsnp > 10) {

                        double localfst = fstdenom > 0? fstnum/fstdenom:0;
                        double divratio = div2>0 ?div1/div2:0;
                        avgcov /= numsnp;
                        ofs << id1 << "\t"<< chr << "\t" << ((long)(window*wnum)) << "\t" << ((long)(window*(wnum+1))) << "\t" << ((long)(window*(wnum+0.5))) << "\t" << numsnp << "\t" << localfst << "\t" << "\t" << div1 << "\t" << div2 << "\t" <<divratio << "\t" << avgcov << endl;
                        fst.push_back (localfst);
                        div.push_back(divratio);
                        cov.push_back (avgcov);
                    }
                    wnum = pos/window;
                    fstnum = 0 ; fstdenom = 0 ;
                    div1 = 0 ; div2 = 0 ;
                    avgcov =  0;
                    numsnp = 0 ;
                }
            }

            if (flag ) { 
                div1 += 2.*a1*(n1-a1)/(n1*(n1-1));
                div2 += 2.*a2*(n2-a2)/(n2*(n2-1));
                double h1 = (1.*a1*(n1-a1))/(n1*(n1-1)); 
                double h2 = (1.*a2*(n2-a2))/(n2*(n2-1)); 
                double tmp = pow( ((1.*a1)/n1-(1.*a2)/n2),2) - h1/n1 - h2/n2;
                if ( f1 >=0.05 && f1<= 0.95 && f2 >=0.05 && f2 <=0.95){
                    fstnum += tmp;
                    fstdenom += (tmp + h1 + h2);
                }
                avgcov += val; 
                numsnp ++;
                y3 = std::getline (cafs,line3).good();
            }

            y1 = std::getline (p1afs,line1).good();
            y2 = std::getline (p2afs,line2).good();
        }
        p1afs.close ();
        p2afs.close ();
        cafs.close();
    }
    ofs.close ();
}

void computexaratio::process ( string pop1, string pop2, string pop3, string covariates, string output, vector<double> &drift1, vector<double> &drift2, vector<double> &drift3, vector<double> &div, vector<double> &cov ){

    vector<string> p1;
    vector<string> p2;
    vector<string> p3;
    vector<string> c;

    functions::tokenize (pop1, p1, ",");
    functions::tokenize (pop2, p2, ",");
    functions::tokenize (pop3, p3, ",");
    functions::tokenize (covariates, c, ",");
    if ( !(p1.size()==p2.size() && p1.size()==c.size())) { 
        cerr << "List of files " << pop1 << "," << pop2  << "," << covariates << " must be equal " << endl;
        exit (1);
    }
    ofstream ofs (output.c_str());

    for (int i  = 0 ;  i < p1.size(); i ++)  {
        if (io::debug>=0)
            cout << "Opening file   " << p1[i] << endl;
        ifstream p1afs (p1[i].c_str());
        ifstream p2afs (p2[i].c_str());
        ifstream cafs (c[i].c_str());
        ifstream p3afs (p3[i].c_str());

        string line1;
        string line2;
        string line3;
        string line4;

        bool y1 = std::getline (p1afs,line1).good();
        bool y2 = std::getline (p2afs,line2).good();
        bool y3 = std::getline (cafs,line3).good();
        bool y4 = std::getline (p3afs, line4).good();

        int wnum  = -1;

        double fstnum12 = 0 ;
        double fstdenom12 = 0 ;
        double fstnum23 = 0 ;
        double fstdenom23 = 0 ;
        double fstnum31 = 0 ;
        double fstdenom31 = 0 ;


        double div1 = 0 ; 
        double div2 = 0 ;
        double avgcov =  0;
        int numsnp = 0 ;
        int linenum = 0 ;
 
        while (y1 && y2 && y3 && y4) { 

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
            while (y3 && line3[0]=='#') { 
                y3 = std::getline (cafs,line1).good();
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
            double genpos;
            if (usegenpos) { 
                if (snpmap->find(id1)!=snpmap->end ()){
                    string l = (*snpmap)[id1];
                    vector<string> ltoks;
                    functions::tokenize (l.c_str(), ltoks, " \t");
                    genpos = atof(ltoks[2].c_str());
                } else {
                    cerr << "Cannot find SNP " << id1 <<endl;
                }        
            }

            toks1.resize(0);
            functions::tokenize (line2.c_str(), toks1, " \t");
            string id2  = toks1[0].c_str();
            int a2 = atoi(toks1[1].c_str()); 
            int b2 = atoi(toks1[2].c_str());
            int n2 = a2 + b2;
            double f2 = (1.*a2)/n2;

            toks1.resize(0);
            functions::tokenize (line3.c_str(), toks1, " \t");
            string id3  = toks1[0];
            double val =  atof(toks1[1].c_str());
            if (io::debug >= 1) 
                cout << "Here " << id1 << "\t" << id3 << endl;
            

            toks1.resize(0);
            functions::tokenize (line4.c_str(), toks1, " \t");
            string id4  = toks1[0].c_str();
            int a3 = atoi(toks1[1].c_str()); 
            int b3 = atoi(toks1[2].c_str());
            int n3 = a3 + b3;
            double f3 = (1.*a3)/n3;


            if ( !(id1.compare(id2)==0 )) {
                cerr << "Files " << p1[i] << "," << p2[i] << "," << c[i] << " differ at line " << linenum << endl;
                exit (1);
            } else  {
                if ( !(id1.compare (id3)==0)){
                    toks2.resize(0);
                    functions::tokenize (id3.c_str(), toks2, ":");
                    long pos3 = atol (toks2[1].c_str());
                    if ( pos3 < pos )  { 
                        y3 = std::getline (cafs,line3).good();
                        continue;
                    } else 
                        flag = false;
                }
            }
            if (io::debug >= 1) 
                cout << "There " << flag  << "\t" << wnum  << endl;

 
            if (wnum == -1 || (wnum + 1)  * window <= pos ) { 
                if (wnum >=0 && numsnp > 10) {

                    double localfst12 = fstdenom12 > 0? fstnum12/fstdenom12:0;
                    double localfst23 = fstdenom23 > 0? fstnum23/fstdenom23:0;
                    double localfst31 = fstdenom31 > 0? fstnum31/fstdenom31:0;
                    double localdrift1 = localfst12 + localfst31 - localfst23;
                    double localdrift2 = localfst12 + localfst23 - localfst31;
                    double localdrift3 = localfst23 + localfst31 - localfst12;
                    double divratio = div2>0 ?div1/div2:0;
                    avgcov /= numsnp;
                    ofs << id1 << "\t"<< chr << "\t" << ((long)(window*(wnum+0.5))) << "\t" << numsnp << "\t" << localdrift1 << "\t" << localdrift2 << "\t" << localdrift3 << "\t" << div1 << "\t" << div2 << "\t" <<divratio << "\t" << avgcov << endl;
                    drift1.push_back (localdrift1);
                    drift2.push_back (localdrift2);
                    drift3.push_back (localdrift3);
                    div.push_back(divratio);
                    cov.push_back (avgcov);
                }
                wnum = pos/window;
                fstnum12 = 0 ; fstdenom12 = 0 ;
                fstnum23 = 0 ; fstdenom23 = 0 ;
                fstnum31 = 0 ; fstdenom31 = 0 ;
                div1 = 0 ; div2 = 0 ;
                avgcov =  0;
                numsnp = 0 ;
            }

            if (flag ) { 
                div1 += 2.*a1*(n1-a1)/(n1*(n1-1));
                div2 += 2.*a2*(n2-a2)/(n2*(n2-1));
                double h1 = (1.*a1*(n1-a1))/(n1*(n1-1)); 
                double h2 = (1.*a2*(n2-a2))/(n2*(n2-1)); 
                double h3 = (1.*a3*(n3-a3))/(n3*(n3-1)); 
                double tmp = pow( ((1.*a1)/n1-(1.*a2)/n2),2) - h1/n1 - h2/n2;
                if ( f1 >=0.05 && f1<= 0.95 && f2 >=0.05 && f2 <=0.95){
                    fstnum12 += tmp;
                    fstdenom12 += (tmp + h1 + h2);
                }


                tmp = pow( ((1.*a3)/n3-(1.*a2)/n2),2) - h3/n3 - h2/n2;
                if ( f3 >=0.05 && f3<= 0.95 && f2 >=0.05 && f2 <=0.95){
                    fstnum23 += tmp;
                    fstdenom23 += (tmp + h3 + h2);
                }

                tmp = pow( ((1.*a3)/n3-(1.*a1)/n1),2) - h3/n3 - h1/n1;
                if ( f3 >=0.05 && f3<= 0.95 && f1 >=0.05 && f1 <=0.95){
                    fstnum31 += tmp;
                    fstdenom31 += (tmp + h1 + h3);
                }

                avgcov += val; 
                numsnp ++;
                y3 = std::getline (cafs,line3).good();
            }

            y1 = std::getline (p1afs,line1).good();
            y2 = std::getline (p2afs,line2).good();
            y4 = std::getline (p3afs,line4).good();
        }
        p1afs.close ();
        p2afs.close ();
        p3afs.close ();
        cafs.close();
    }
    ofs.close ();
}


int main (int argc, char *argv[]) { 
    computexaratio c (argc, argv);
}
