#include "ram.h"
#include "snp.h"
#include "data.h"
#include <boost/heap/fibonacci_heap.hpp>
using namespace boost::heap;

struct cmp  { 
    bool operator () (const pair<string, double> &a, const pair<string, double> &b) const {
        return a.second < b.second;
    }
};


class ranker {
    public:
        ranker (int);
        ranker (int argc, char *argv[]) ;
        void read_freq ();
        void pick_tagged  (string id) ;
        void read_r2map () ;
        double get_rank (string id) ;

        data *d ,*p;
        string rankfile;
        string freqfile;
        string featurefile;
        string r2file;
        string confidencefile;
        double r2threshold;
        double freqthreshold;
        bool onlyconfident;
        double confidentweight;
        bool givenotherinfo;
        string otherinfofile;

        unordered_map <string, int> *confidencemap;
        unordered_map <string, int> *featuremap;
        unordered_map <string, double> *freqmap; 
        unordered_map <string, double> *r2map;
        unordered_map < string ,vector < string> > *tagmap;
        unordered_map < string, bool> *ignoremap;
        unordered_map < string, bool> *donemap;
        unordered_map < string, double> *tagsizemap;
        unordered_map < string, double> *rankmap;
        unordered_map < string, double> *yrimap;

        int rankingmethod;
        int nsnp; 
        double nind;

        fibonacci_heap < pair<string, double> , boost::heap::compare<cmp> >  heap;
        typedef fibonacci_heap< pair<string, double> , boost::heap::compare<cmp> > ::handle_type handle_h;
        unordered_map <string, handle_h> *handlemap;
};


ranker::ranker (int x)  {
    fibonacci_heap<int> fib_heap;

    fib_heap.push(2);
    fib_heap.push(1);
    fib_heap.push(3);
    cout << "Here" << endl;
    for (fibonacci_heap<int>::ordered_iterator i = fib_heap.ordered_begin (); i != fib_heap.ordered_end (); i++) { 
        cout << *i  << endl;
    }

    fibonacci_heap < pair<string, double> , boost::heap::compare<cmp> >  h;
    handlemap = new unordered_map <string, handle_h > ();
    (*handlemap)["a"]=h.push( pair<string, double> ("a",1));
    (*handlemap)["b"]=h.push( pair<string, double> ("b",2));
    (*handlemap)["c"]=h.push( pair<string, double> ("c",3));
    (*handlemap)["d"]=h.push( pair<string, double> ("d",4));
    (*handlemap)["e"]=h.push( pair<string, double> ("e",5));

    for (fibonacci_heap<pair<string, double> , boost::heap::compare<cmp> >::ordered_iterator i = h.ordered_begin (); i != h.ordered_end (); i++) { 
        cout << i->first << "\t" << i->second  << endl;
    }

    cout << "Here" << endl;
    int j  = 0 ;
    while (h.size () >  0) {
        pair<string, double> p = h.top ();
        h.pop ();
        cout << p.first << "\t" << p.second << endl;
        j ++;
        if (j==1) {
            handle_h &ht1= (*handlemap)["b"]; (*ht1).second = 4;
            handle_h &ht2= (*handlemap)["c"]; (*ht2).second = 0;
            h.update(ht1);
            h.update(ht2);
        }
    
    }

    /*
    for (fibonacci_heap<pair<string, double> , boost::heap::compare<cmp> >::ordered_iterator i = h.ordered_begin (); i != h.ordered_end (); i++) { 
        cout << i->first << "\t" << i->second  << endl;
    }*/
}

ranker::ranker (int argc, char *argv[]) { 
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

    freqfile =  ""; r2file = "" ; featurefile  = ""; confidencefile = "";
    p->get_string ("rankfile", rankfile, true);
    p->get_string ("freqfile", freqfile, true);
    p->get_string ("r2file", r2file, true);
    p->get_string ("featurefile", featurefile, false);
    p->get_string ("confidencefile", confidencefile, false);
    p->get_double ("r2threshold", r2threshold, 0.8, false);
    p->get_double ("freqthreshold", freqthreshold, -1, false);
    p->get_boolean ("onlyconfident", onlyconfident,  false, false);
    p->get_double ("confidentweight", confidentweight, 10, false);
    p->get_int ("debug", io::debug, 0);
    p->get_int ("rankingmethod", rankingmethod, 0, false);
    givenotherinfo = p->get_string ("otherinfo", otherinfofile, false);
    nind = -1;

    rankmap = new unordered_map <string, double> ();
    handlemap = new unordered_map <string, handle_h > ();
    featuremap = new unordered_map <string, int> ();
    yrimap = new unordered_map <string, double> ();
    fileio::read_map (otherinfofile, yrimap, 0, 1);

    fileio::read_map (featurefile, featuremap, 0 ,1);

    confidencemap = new unordered_map <string, int> ();
    fileio::read_map (confidencefile, confidencemap, 0 ,1);

    read_freq ();

    if (onlyconfident) { 
        for (unordered_map<string,int>::iterator i = featuremap->begin(); i != featuremap->end(); i++) { 
            if ((confidencemap->find(i->first)==confidencemap->end()) || ((*confidencemap)[i->first]==0))
                (*ignoremap)[i->first] = true;
        }
    }

    if (givenotherinfo) {
        for (unordered_map <string,double>::iterator i = freqmap->begin(); i != freqmap->end(); i++) { 
            if ( i->second < 3*(*yrimap)[i->first] ) 
                (*ignoremap)[i->first] = true;
        }
    }
    read_r2map ();

// Read all inputs

    for (unordered_map<string,int>::iterator i = featuremap->begin(); i != featuremap->end(); i++) { 
        string id = i->first;
        vector<string> &tags = (*tagmap)[id];
        int d = 1 ;
        for (int j = 0  ; j < tags.size(); j++) { 
            if (!(*ignoremap)[tags[j]])
                d ++;
        }
        (*tagsizemap)[id] = d;
    }

    if (io::debug >= 1){ 
        for (unordered_map<string,int>::iterator i = featuremap->begin(); i != featuremap->end(); i++) { 
            cout << i->first << "\t" << (*ignoremap)[i->first] << "\t" << (*confidencemap)[i->first] << "\t" << (*tagsizemap)[i->first] << endl;
        }
    }
    typedef fibonacci_heap<pair<string, double> , boost::heap::compare<cmp> > fheap;

    for (unordered_map<string,int>::iterator i = featuremap->begin(); i != featuremap->end(); i++) { 
        string id = i->first;
        double rank = get_rank (id);
        (*handlemap)[id]=heap.push (pair<string, double> (id,rank));
    }

    ofstream rankfs ( rankfile);
    int count  = 0 ;
    while (heap.size()>1) {
        count++;
        pair<string, double> p = heap.top ();  heap.pop ();
        string id = p.first;
        double rank  = p.second;
        if (io::debug>=0)
            cout << id << "\t" << rank << "\t" << heap.size() <<endl;

        if ((*donemap)[id])
            continue;
        (*donemap)[id] = true;
        if ((*ignoremap)[id])
            continue;

        vector<string> &tags = (*tagmap)[id];
        rankfs << id << "\t" << (*confidencemap)[id] << "\t" << (*freqmap)[id] << "\t" << (*featuremap)[id] << "\t" << (*tagsizemap)[id] << "\t" << rank << endl;
        unordered_set <string> changeset;

        for (int j = 0  ; j < tags.size(); j++) { 
            string id1 = tags[j];
            (*donemap)[id1] = true;
            vector<string> &tags1 = (*tagmap)[id1];
            for (int k  =0  ; k < tags1.size(); k++){
                string id2 = tags1[k];
                if ((*donemap)[id2]) continue;
                if ((*ignoremap)[id2]) continue;
                if (id2.compare(id)!=0) {
                    (*tagsizemap)[id2]=(*tagsizemap)[id2]-1;
                    changeset.insert(id2);
                }
            }
        }

        /*
        handle_h &ht = (*handlemap)[id];
        (*ht).second = -1;
        heap.update (ht);*/
        if (io::debug >= 0) {
            cout << "changeset (" << changeset.size () << ") = "   << flush;
            for (unordered_set<string>::iterator itr = changeset.begin (); itr != changeset.end(); itr++) { 
                string id1 = *itr;
                cout <<  id1 << "\t";
            }
            cout << endl;
        }
        for (unordered_set<string>::iterator itr = changeset.begin (); itr != changeset.end(); itr++) { 
            string id1 = *itr;
            double rank = get_rank (id1);
            handle_h &ht1 = (*handlemap)[id1];
            (*ht1).second = rank;
            heap.update (ht1);
        }
    }
}


double ranker::get_rank (string id) { 
    double v = (*tagsizemap)[id];
    double f = (*featuremap)[id];
    double p = (*freqmap)[id];
    int confident = 0;
    double rank =  0;
    if (confidencemap->find(id)!= confidencemap->end())
        confident = (*confidencemap)[id];
    switch(rankingmethod) { 
        case 0:
            rank = v;
            break;
        case 1:
            rank  = p*v/f;
            break;
        case 2:
            rank  = p*v/f;
            rank  = confident * confidentweight   + rank;
            break;
    }
    return rank;
}

void ranker::pick_tagged  (string id) {
    vector<string>& tags = (*tagmap) [id];
    for (int i = 0 ;  i < tags.size(); i++) { 
        string tmpid = tags[i];
        vector<string> &tmptags = (*tagmap)[tmpid];

    }
}

void ranker::read_freq (){ 
	ifstream inp (freqfile.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< freqfile <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
    freqmap = new unordered_map <string, double> ();
    ignoremap = new unordered_map <string, bool> ();
    tagsizemap = new unordered_map <string, double> ();
    donemap = new unordered_map <string, bool> ();
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
        vector<string> toks;
        functions::tokenize(line.c_str(),toks," \t");
        string val = line;
        nind = atof(toks[1].c_str()) + atof(toks[2].c_str());
        double freq =(1.* atof(toks[1].c_str()) )/nind;
        string id = toks[0];
        (*freqmap)[id] = freq;
        if (freq < freqthreshold ) 
            (*ignoremap)[id] = true;
        else
            (*ignoremap)[id] = false;
        if ( io::debug>=1) 
            cout << id << "\t" << freq << "\t" << freqthreshold << "\t" << (*ignoremap)[id] << endl;
        (*tagsizemap)[id] = 1 ; 
        (*donemap)[id] = false;
        nsnp++;
	}	
	inp.close ();
}

void ranker::read_r2map () { 
   	ifstream inp (r2file.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< r2file <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
    r2map = new unordered_map < string , double >  ();
    tagmap = new unordered_map <string, vector<string> > ();
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
        vector<string> toks;
        functions::tokenize(line.c_str(),toks," \t");
        string id1 = toks[0] + ":" + toks[1]; 
        string id2 = toks[2] + ":" + toks[3];
        if (onlyconfident) {
            if (confidencemap->find (id1) == confidencemap->end ())
                continue;
            if (confidencemap->find (id2) == confidencemap->end ())
                continue;
            int c1 = (*confidencemap)[id1];
            int c2 = (*confidencemap)[id2];
            if (c1==0 || c2==0)
                continue;
        }
        double r2 = atof(toks[toks.size()-1].c_str());
        if (r2 < r2threshold ) 
            continue;
        if ((freqmap->find (id1) == freqmap->end())||(freqmap->find(id2) == freqmap->end()))
            continue;

        string  k = id1 + "," + id2;
        (*r2map)[k] = r2;
        if (tagmap->find(id1)==tagmap->end ()){ 
            (*tagmap)[id1] = vector<string> ();
        }
        vector<string> &tmp1 = (*tagmap)[id1];
        tmp1.push_back (id2);

        if (tagmap->find(id2)==tagmap->end ()){ 
            (*tagmap)[id2] = vector<string> ();
        }
        vector<string> &tmp2 = (*tagmap)[id2];
        tmp2.push_back (id1);

        (*tagsizemap)[id1] = (*tagsizemap)[id1]+1;
        (*tagsizemap)[id2] = (*tagsizemap)[id2]+1;

	}	
	inp.close ();
}

int main (int argc, char *argv[]) {
    ranker r (argc,argv);
//    ranker r (0);
}
