#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"

class polarize { 
	public:
		polarize (int argc, char *argv[]) ;
		void convert ();

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

};


ofstream outfs;

polarize::polarize (int argc, char *argv[]) {
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
	
}

void polarize::convert ()  {

	string polarizefile;
	string freqoutfile;
	bool polarize = p->get_string ("polarizesnpname", polarizefile, false);
	p->get_string("freqoutname|freqoutfilename",freqoutfile,true);
	ofstream freqfs (freqoutfile.c_str());

	unordered_map<string, string > *pmapsnps = NULL;
	if (polarize){  
		pmapsnps = new unordered_map <string, string> ();
		fileio::read_map (polarizefile,  pmapsnps, 0 , 1);
		if (io::debug >= 2)  {
			typedef unordered_map<string, string>::iterator iter;
			for ( iter i  = pmapsnps->begin();  i != pmapsnps->end(); i++){
				cout << i->first << "\t" << i->second << endl;
			}
		}
	}
	string freqfile;
	p->get_string ("freqname", freqfile, true);

	
	ifstream inp (freqfile.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< freqfile <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
	
		istringstream ss (line);

		string id;
		double vfreq;
		double rfreq;
		string var;
		string ref;

		ss >> id;
		ss >> vfreq;
		ss >> rfreq;
		ss >> var;
		ss >> ref;


		if (polarize) {
				if ( pmapsnps->find (id) != pmapsnps->end ()) {
					string a = (*pmapsnps)[id];
					if (io::debug>= 2) {
						cout << "Reference " << id << "\t" << a << endl;
					}

					if (a.compare("?")==0)
						continue;
					else { 
						if (ref.compare("X")!=0 && var.compare("X")!=0) { 
							if (var.compare(a)==0){
								var = ref;
								ref = a;
								double tmp = vfreq;
								vfreq = rfreq;
								rfreq = tmp; 
							}
						} else  {
							if (ref.compare("X")==0) {
								if (var.compare(a)==0){
									ref = a;
									var = "X";
									double tmp = vfreq;
									vfreq = rfreq;
									rfreq = tmp; 
								} else {
									ref = a;
								}
							} else if (var.compare ("X")==0){
								if (ref.compare (a)==0) {
								} else {
									var = ref;
									ref = a;
									double tmp = vfreq;
									vfreq = rfreq;
									rfreq = tmp; 
								}
							}
						}
					}

				} else 
					continue;
			}
		freqfs << id << "\t" << vfreq << "\t" << rfreq << "\t" << var <<"\t" <<ref <<endl;
	}
	inp.close ();
	freqfs.close ();

}


int main (int argc, char* argv[]) {
	polarize c = polarize(argc,argv);
	c.convert ();
}
