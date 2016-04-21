#include "ram.h"


class annotate { 

    public:
        annotate  (int argc, char *argv[]) ;
    vector<int> start; vector<int> end; vector<int> val;
    vector<string> sval;
    data *d;
    data *p;
    bool stringval;
};

annotate::annotate  (int argc, char *argv[]) { 
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
 
    string f1;
    p->get_string ("snpfile",f1,true);
    string f2;
    p->get_string ("intervalfile", f2, true);
    string f3;
    p->get_string ("annotatedfile", f3, true);
    p->get_boolean ("stringval", stringval, false, false);

    vector<string> *tmp1 = new vector<string> ();
    fileio::read_vector ( f2, tmp1, -1 );
    for (int i  = 0 ;  i < tmp1->size(); i++) { 
        string s = (*tmp1)[i];
        vector<string> toks;
        functions::tokenize (s, toks, " \t");
        start.push_back (atoi(toks[1].c_str()));
        end.push_back (atoi(toks[2].c_str()));
        val.push_back (atoi(toks[3].c_str()));
        if (stringval)
            sval.push_back (toks[3]);

    }


	ifstream inp (f1.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< f1 <<endl;
		exit(1);
	}

    string line;
	int linenum  = 0;
    int idx =  0;
    int limit = start.size();

    ofstream ofs ( f3.c_str());
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;

		istringstream ss (line);
     	string id;
		string chr;
		double genpos;
		double physpos;
		string var;
		string ref;

		ss >> id;
		ss >> chr;
		ss >> genpos;
		ss >> physpos;
		ss >> var;
		ss >> ref;

        if (io::debug>=1)
            cout << physpos << "\t" << idx << "\t" << start[idx] << "\t" << end[idx] << endl;
        if (physpos  < start[idx]) {
        } else {
           while (idx <  limit && physpos > end[idx]) idx++;
           if ( idx < limit)  { 
                if (stringval)
                    ofs << id << "\t" << sval[idx] << endl; 
                else
                    ofs << id << "\t" << val[idx] << endl; 
           }
        }
    }
    inp.close ();
    ofs.close ();
	
}


int main (int argc, char *argv[]) { 
    annotate a (argc, argv);
}

