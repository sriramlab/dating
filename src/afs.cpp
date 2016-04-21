#include "ram.h"
#include "afs.h"


afs::afs ()  {
    alloced =false;
    npops  = 0 ;
    counts  = NULL;
    countsize = 1;
}

afs::afs (string filename)  {
    alloced =false;
    npops  = 0 ;
    counts  = NULL;
    countsize = 1;
    read(filename);
}

afs::afs (vector<int> &p, vector<string> &pl) { 
    npops = p.size();
    if ( npops != pl.size()){
        cerr << "Invalid population sizes " << endl;
        exit (1);
    }
    chrsperpop.resize (npops);
    poplabels.resize (npops);
    countsize = 1;
    for (int i  = 0 ; i < npops; i++) { 
        chrsperpop[i] = p[i];
        poplabels[i] = pl[i];
        popmap[poplabels[i] ] = i;
        countsize *= (chrsperpop[i]+1);
    }
    counts = new double [countsize];
    for (int i =  0 ;  i < countsize; i++)
        counts[i] =  0;
    alloced = true;
}

afs* afs::permute (string labels){ 
    size_t found = labels.find ("sampsizes:");
    if (found != std::string::npos) {
        labels.erase (found, found + 10);
    } 
    boost::erase_all (labels, " ");

    vector<string> toks;
    vector<string> toks1;

    functions::tokenize (labels, toks, ";");
    vector<string> pops (toks.size());
    for (int  i = 0 ; i < toks.size(); i++){ 
        toks1.resize(0);
        functions::tokenize (toks[i], toks1,":");
        pops[i] = toks1[0];
    }
    for (int i = 0  ; i < pops.size(); i++)
        cout << pops[i] << "\t";
    cout << endl;

    return permute (pops);
}

afs* afs::permute (vector<string> labels) { 

    if (labels.size()!= npops){
        cerr << "Error in permute :" <<  labels.size() << "!= " << npops <<endl;
        exit(1);
    }
    // perm
    vector<int>  perm;
    for (int i = 0 ; i < labels.size();i++){
        if (popmap.find(labels[i])==popmap.end()){
            cerr << "Error in permute : cannot find " << labels[i] << endl;
            exit(1);
        }
        perm.push_back (popmap[labels[i]]);
    }
    vector<int> tmpsizes (npops);
    for (int i = 0 ; i < npops; i++){  
        tmpsizes[i] = chrsperpop[perm[i]];
    }

    afs *a = new afs (tmpsizes, labels);

    vector<int> counter (npops);
    vector<int> tmpcounter (npops);
    while (1) { 
        /*
        string s = "";
        for (int i  = 0 ; i < counter.size(); i++)
            s += tostring(counter[i]) + "\t";
        s += tostring (get(counter) );
        cout << s << endl;*/

        for (int j = 0 ; j < npops ; j++)
            tmpcounter [ perm[j] ] = counter[j];

        a->put (counter, get(tmpcounter));

        int i = npops - 1;
        for (; i>=0 && counter[i] == tmpsizes[i] ; i--) ;
        if (i<0)
            break;
        counter[i]++;
        for (int j  = i+1; j < counter.size(); j++)
            counter[j] =  0;

    }
    return a;
}

afs::afs (vector<int> &p) { 
    npops = p.size();
    chrsperpop.resize (npops);
    poplabels.resize (npops);
    countsize = 1;
    for (int i  = 0 ; i < npops; i++) { 
        chrsperpop[i] = p[i];
        poplabels[i] = tostring(i);
        popmap[poplabels[i] ] = i;
        countsize *= (chrsperpop[i]+1);
    }
    counts = new double [countsize];
    for (int i =  0 ;  i < countsize; i++)
        counts[i] =  0;

    alloced = true;
}

afs::~afs () {
    if (alloced ) {
        chrsperpop.clear ();
        poplabels.clear ();
        delete[] counts;
    }
}

inline double afs::get (vector<int> &index) { 
    int tmp  = index[0] ;
    for (int i  = 1 ; i < npops; i++) { 
        tmp *= (chrsperpop[i]+1);
        tmp += index[i];
    }
    return counts[tmp];
}

inline void afs::put (vector<int> &index, double v) { 
    int tmp  = index[0] ;
    for (int i  = 1 ; i < npops; i++) { 
        tmp *= (chrsperpop[i]+1);
        tmp += index[i];
    }
    counts[tmp] = v;
}

void afs::increment (vector<int> &index) { 
    int tmp  = index[0] ;
    for (int i  = 1 ; i < npops; i++) { 
        tmp *= (chrsperpop[i]+1);
        tmp += index[i];
    }
    counts[tmp] =  counts[tmp]  + 1;

}

void afs::increment (vector<int> &index, double v) { 
    int tmp  = index[0] ;
    for (int i  = 1 ; i < npops; i++) { 
        tmp *= (chrsperpop[i]+1);
        tmp += index[i];
    }
    counts[tmp] =  counts[tmp]  + v;

}

void afs::increment (int *index, int n) { 
    if ( n != npops)
        exit(1);
    vector<int> tmpindex (npops);
    for (int i  = 0 ; i < npops; i++)
        tmpindex[i] =  index[i];
 
    increment (tmpindex);
}


inline double afs::operator () (vector<int> &index ) {
    return get (index);
}

inline double afs::operator () (int *index, int n ) {
    if ( n != npops)
        exit(1);
    vector<int> tmpindex (npops);
    for (int i  = 0 ; i < npops; i++)
        tmpindex[i] =  index[i];
    return get (tmpindex);
}




void afs::read (string filename ) {  
    ifstream ifs (filename.c_str());
    string line ;
    int count =  0;
    bool header = true;
    while ( std::getline (ifs,line) ) { 
        char c = line[0];
        if (header && c=='#') { 
            if (alloced) { 
                count ++;
                if (count == 1){
                    vector<string> toks;
                    functions::tokenize (line, toks, " \t");
                    if (toks.size()!=npops+1) {
                        cerr << "Invalid header in " << filename << endl;
                        exit(1);
                    }
                    for (int i  = 1; i < toks.size(); i++) {
                        poplabels[i-1] = toks[i];
                        popmap[poplabels[i-1] ] = i-1;
                    }
                }
            } else { 
                count ++;
                vector<string> toks;
                functions::tokenize (line, toks, " \t");
                if (count == 1){
                    npops = toks.size()-1;
                    chrsperpop.resize (npops);
                    poplabels.resize (npops);
                    for (int i  = 1; i < toks.size(); i++) {
                        poplabels[i-1] = toks[i];
                        popmap[poplabels[i-1] ] = i-1;
                    }
                } else if (count==2){
                    for (int i  = 1; i < toks.size(); i++) {
                        chrsperpop[i-1] = atoi(toks[i].c_str());
                        countsize *= (chrsperpop[i-1]+1);
                    }
                    counts = new double [countsize];
                    for (int i =  0 ;  i < countsize; i++)
                        counts[i] =  0;
                    alloced = true;
                } 
            }
        } else {
            header =false;
            if (alloced ) { 
                vector<string> toks;
                functions::tokenize (line, toks, " \t");
                if (toks.size ()!=npops+1) { 
                    cerr << "Invalid line in " << filename << endl;
                    exit (1);
                }
                vector<int> index (npops);
                for (int i = 0 ;  i < npops ; i++)
                    index[i] = atoi(toks[i].c_str());
                put (index,atof(toks[toks.size()-1].c_str()));

            } else { 
            }
        }
    }
    ifs.close ();
}

afs* afs::ascertain (vector<string> &pop, vector<int> &dc, vector<int>  &tc) { 
    vector<int> newsizes ;
    int newnpops;
    vector<string> newlabels;
    // -1 if there is no ascertianment in a population
    vector<int> ascertaindc(npops,-1);
    vector<int> ascertainac(npops,-1);
    for (int i = 0 ;  i < npops; i++) {
        int k  = -1;
        for (int j =  0 ; j < pop.size(); j++){ 
            if (poplabels[i].compare(pop[j])==0) {
                k = j;
                break;
            }
                
        }
        if ( k>=0) { 
            // ascertainment in population i
            if ( tc[k] > chrsperpop[i]) {
                cerr << "Invalid ascertainment" << endl;
                exit(1);
            } else if (tc[k] == chrsperpop[i]) {
            } else {
                int n = chrsperpop[i] - tc[k];
                newsizes.push_back ( n );
                newlabels.push_back (poplabels[i]);
            }
            ascertaindc[i] = dc[k];
            ascertainac[i] = tc[k]-dc[k];
        } else { 
            // no ascertainment in pop i
            newsizes.push_back (chrsperpop[i]);
            newlabels.push_back (poplabels[i]);

        }
    }

    newnpops = newsizes.size();
    afs *a = new afs (newsizes, newlabels);

    vector<int> counter (npops);
    vector<int> acounter (npops);
    vector<int> tmpcounter(newnpops);
    while (1) { 
        /*
        string s = "";
        for (int i  = 0 ; i < counter.size(); i++)
            s += tostring(counter[i]) + "\t";
        s += tostring (get(counter) );
        cout <<  "s =  "  <<  s << endl;
        */

        bool flag = true;
        tmpcounter.resize(newnpops);
        double lnprob = 0;
        for (int i = 0,k=0 ; i < counter.size(); i++)  {
            acounter[i] = chrsperpop[i] - counter[i];

            // Ascertain pop i
            if (ascertaindc[i]>=0){ 
                // Cannot ascertain this configuration
                if (counter[i] < ascertaindc[i] || acounter[i] < ascertainac[i]) {
                    flag = false;
                    break;
                } else {
                    // After ascertainment, there are no more alleles left. Drop
                    // this population
                    if (counter[i]==ascertaindc[i] && acounter[i] == ascertainac[i])
                        continue;
                    double tmpp= mathfn::lnchoose(counter[i],ascertaindc[i]) +mathfn::lnchoose(acounter[i],ascertainac[i])-mathfn::lnchoose(chrsperpop[i],(ascertaindc[i]+ascertainac[i]));
                    lnprob += tmpp;
                    // If there are alleles left,include this population
                    tmpcounter[k] = counter[i] - ascertaindc[i];
                    k++;
                }
            } else {
                tmpcounter[k] = counter[i];
                k++;
            }
        }
        
        if (flag ) {
            /*
            s = "";
            for (int i  = 0 ; i < tmpcounter.size(); i++)
                s += tostring(tmpcounter[i]) + "\t";
            cout << "tmp counter = " << s << "\t" << get(counter) << endl;
            */

            a->increment (tmpcounter, get(counter)*exp(lnprob));
        }

        int i = npops - 1;
        for (; i>=0 && counter[i] == chrsperpop[i] ; i--) ;
        if (i<0)
            break;
        counter[i]++;
        for (int j  = i+1; j < counter.size(); j++)
            counter[j] =  0;

    }
    return a;
}

afs* afs::ascertain (string asc) { 
    size_t found = asc.find ("ascertain:");
    if (found != std::string::npos) {
        asc.erase (found, found + 10);
    } 
    boost::erase_all (asc, " ");

    vector<string> toks;
    vector<string> toks1;

    functions::tokenize (asc, toks, ";");
    vector<string> pops (toks.size());
    vector<int> dc (toks.size());
    vector<int> tc (toks.size());
    for (int  i = 0 ; i < toks.size(); i++){ 
        toks1.resize(0);
        functions::tokenize (toks[i], toks1,"::");
        pops[i] = toks1[0];
        dc[i] = atoi(toks1[1].c_str());
        tc[i] = atoi(toks1[2].c_str());
    }

    if (io::debug >= 1) { 
        for (int i  = 0 ; i < pops.size(); i++) { 
            cout << pops[i]  << "," << dc[i] << "\t" << tc[i] << endl;
        }
    }
    return ascertain (pops, dc, tc);
}



/*
afs* afs::subsample (vector<string> &pop, vector<int>  &tc) { 
    vector<int> newsizes ;
    int newnpops;
    vector<string> newlabels;
    // -1 if there is no ascertianment in a population
    vector<int> ascertaindc(npops,-1);
    vector<int> ascertainac(npops,-1);
    for (int i = 0 ;  i < npops; i++) {
        int k  = -1;
        for (int j =  0 ; j < pop.size(); j++){ 
            if (poplabels[i].compare(pop[j])==0) {
                k = j;
                break;
            }
                
        }
        if ( k>=0) { 
            // ascertainment in population i
            if ( tc[k] > chrsperpop[i]) {
                cerr << "Invalid ascertainment" << endl;
                exit(1);
            } else if (tc[k] == chrsperpop[i]) {
            } else {
                int n = chrsperpop[i] - tc[k];
                newsizes.push_back ( n );
                newlabels.push_back (poplabels[i]);
            }
            ascertaindc[i] = dc[k];
            ascertainac[i] = tc[k]-dc[k];
        } else { 
            // no ascertainment in pop i
            newsizes.push_back (chrsperpop[i]);
            newlabels.push_back (poplabels[i]);

        }
    }

    newnpops = newsizes.size();
    afs *a = new afs (newsizes, newlabels);



} */

void afs::print (string filename)  {
    ofstream ofs (filename.c_str());
    ofs << "#\t" ;
    std::copy (poplabels.begin(), poplabels.end(), ostream_iterator<string>(ofs, "\t"));
    ofs << endl;
    ofs << "#";
    for (int i = 0 ; i < chrsperpop.size(); i++)
        ofs << "\t" << chrsperpop[i] ;
    ofs << endl;

    vector<int> counter (npops);
    while (1) { 
        string s = "";
        for (int i  = 0 ; i < counter.size(); i++)
            s += tostring(counter[i]) + "\t";
        s += tostring (get(counter) );
        ofs << s << endl;


        int i = npops - 1;
        for (; i>=0 && counter[i] == chrsperpop[i] ; i--) ;
        if (i<0)
            break;
        counter[i]++;
        for (int j  = i+1; j < counter.size(); j++)
            counter[j] =  0;

    }
    ofs.close ();
}

void afs::print ()  {
    cout << "#\t" ;
    std::copy (poplabels.begin(), poplabels.end(), ostream_iterator<string>(cout, "\t"));
    cout << endl;

    vector<int> counter (npops);
    while (1) { 
        string s = "";
        for (int i  = 0 ; i < counter.size(); i++)
            s += tostring(counter[i]) + "\t";
        s += tostring (get(counter) );
        cout << s << endl;


        int i = npops - 1;
        for (; i>=0 && counter[i] == chrsperpop[i] ; i--) ;
        if (i<0)
            break;
        counter[i]++;
        for (int j  = i+1; j < counter.size(); j++)
            counter[j] =  0;

    }
}



/*
int main (int argc, char *argv[]) { 
    vector<int> popsize (3);
    popsize[0] = 1; popsize[1]=  1; popsize[2] = 3;
    afs a ( popsize);
    a.read ("tmp");
    
    vector<int> tmp (3); tmp[0] = 0; tmp[1] = 1; tmp[2] = 3;
    a.increment (tmp);
    a.print ();
 
    vector<string> ascpops(3); ascpops [0] =  "n"; ascpops[1] = "d"; ascpops [2] = "y";
//    vector<string> ascpops(2); ascpops [0] =  "d"; ascpops[1] = "n";
    vector<int> ad(3); ad[0] =  0; ad[1] = 1    ; ad[2] = 1;
    vector<int> aa(3); aa[0] =  1; aa[1] = 0    ;aa[2] = 1;
    afs *b = a.ascertain  ( ascpops, ad , aa);
    b->print ();
} */  
