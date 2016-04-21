#include "std.h"
#include "data.h"
#include "io.h"
#include "functions.h"

data::data (string filename){
	char **argv = new char *[1];
	argv[0] = new char[filename.length()+1];
	strcpy (argv[0],filename.c_str());
	read_config(1,argv,0);
}

data::data (string filename, map<string, string> &cmap){
	char **argv = new char *[1];
	argv[0] = new char[filename.length()+1];
	strcpy (argv[0],filename.c_str());
	read_config(1,argv,0);
    for (map<string,string>::iterator i = cmap.begin(); i!=cmap.end(); i++){
        configMap[i->first] = i->second;
    }

}

data::data (int argc, char *argv[], const char *optString, const option *longOpts) { 
    if (optString == NULL || longOpts==NULL) {
    optString = "vh";
    const struct option tmplongOpts[] = {
        { "parameter", required_argument, NULL, 'p' },
        { "verbose", no_argument, NULL, 'v' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };
    longOpts = tmplongOpts;


    }

	int opt = 0;
	int longIndex ;
	opt = getopt_long_only( argc, argv, optString, longOpts, &longIndex );
	while( opt != -1 ) {
		if (io::debug >= 4) {
			cout << "opt = " << (char)opt << "\t optind = " << optind << "\t" << optarg << "\t" << longIndex << endl;
		}
		switch( opt ) { 
			case 'v':
				io::debug  = 4;
			break;

			
			case 'h':  
				cout << "Usage:"<<endl;
				exit(1);
				break;
				
			default:
                if (longOpts[longIndex].has_arg)
    				clMap[string(longOpts[longIndex].name)] = string(optarg);
                else
    				clMap[string(longOpts[longIndex].name)] = string("1");
				break;
		}

		opt = getopt_long_only( argc, argv, optString, longOpts, &longIndex );
	}

	read_config (argc,argv,optind);
	
}

void data::read_config (int argc, char *argv[], int optind){
	if (optind < argc) {
		ifstream inp (argv[optind]);
		if (!inp.is_open()){
			cerr << "Error reading file "<< argv[optind] <<endl;
			exit(1);
		}
		string line;
		while ( std::getline (inp, line)){
			char c = line[0];
			if (c=='#')
				continue;
			istringstream ss (line);
			string key;
			string value;
			ss >> key;
			ss >> value;
			char k = key.at(key.length()-1);
			if (k==':'){
				key.erase (key.length()-1,1);
			}
			configMap[key] = value;
			io::println ("key:" + key + "=" + value, 4);
		}
		inp.close();
	}

	map<string,string>::iterator i;
	for ( i = clMap.begin (); i!=clMap.end (); i++) {
		string k = i->first;
		string v = i->second;
		configMap[k] = v;
	}
	macro_sub ();
	
	if (io::debug >=4)
			print_parameters();

}

void data::print_parameters ()  {
	io::println ("###########Parameters#############",0);
	map<string,string>::iterator i;
	for (i = configMap.begin(); i!=configMap.end(); i++){
		string s = "#" + i->first +"\t" + i->second ;
		io::println (s,0);
	}
	io::println ("#################################",0);

}

void data::macro_sub ()  {
	vector<string> macrok;
	vector<string> macrov;
	map <string, int> mmap;	
	map<string,string>::iterator i;
	for (i = configMap.begin(); i!=configMap.end(); i++){
		string k  = i->first;
		bool flag = true;
		for (int j = 0 ; j < k.length(); j++){
			if (islower (k[j]))
				flag = false;
		}
		if (flag) {
			macrok.push_back (k);
			macrov.push_back(i->second);
			mmap[k]=1;
		}
	}

	for (int i = 0; i < macrok.size(); i++) {
		io::println (macrok[i] + "\t" + macrov[i],4);
	}

	for (i = configMap.begin(); i!=configMap.end(); i++){
		string k  = i->first;
		string v = i->second;
		if (mmap.find (k) != mmap.end())
			continue;
		while (1) {
			bool flag = false;
			for (int i = 0 ; i < macrok.size (); i++){
				int f = v.find (macrok[i]);
				if (f!=string::npos) { 
					flag = true;
					v.replace(f,macrok[i].length(),macrov[i]);
				}
			}
			if (!flag)
				break;
		}
		configMap[k] = v;
		bool flag = true;
	}
}





void data::config_simple () {

	io::log = 0;
	if (configMap.find ("logfile") != configMap.end ()) {
		io::logfs.open(configMap["logfile"].c_str());
		io::log = 1;
	}


	
	get_string ("snpname", snpfile, true);
	get_string ("genotypename", genofile, true);
	get_string ("indivname", indfile, true);
	get_double ("binsize", binsize, 0.00001);
	get_string ("output", outputname,"out");
	get_double ("maxdis", maxdis, 1);

	get_int ("debug", io::debug, 0);
	get_int ("seed", seed, 0);
}


bool data::get_boolean (const char *k, bool &var, bool init, bool required) {
	string key(k);
	if (configMap.find (key) != configMap.end()){
		int tmp  = atoi(configMap[key].c_str());
		var = (tmp != 0 );
		return true;
	} else  {
		if (!required) {
			var = init;
			return false;
		} else {
			cerr << key << " required" << endl;
			exit(1);
		}
	}
}

bool data::get_int (const char *k, int &var, int init, bool required) {
	string key(k);
	if (configMap.find (key) != configMap.end()){
		var = atoi(configMap[key].c_str());
		return true;
	} else  {
		if (!required) {
			var = init;
			return false;
		} else {
			cerr << key << " required" << endl;
			exit(1);
		}
	}
}

bool data::get_double_vector (const char *k, vector<double> &v, bool required) {
	string key(k);

	if (configMap.find (key) != configMap.end()){
		vector<string> toks;
		functions::tokenize(configMap[key],toks,",");
		for(int i= 0 ; i < toks.size();i++)
			v.push_back(atof(toks[i].c_str()));
		return true;	
	} else  {
		if (!required) {
			return false;
		} else {
			cerr << key << " required" << endl;
			exit(1);
		}
	}
}

bool data::get_double (const char *k, double &v, double init, bool required) {
	string key(k);

	if (configMap.find (key) != configMap.end()){
		v  = atof(configMap[key].c_str());
		return true;	
	} else  {
		if (!required) {
			v = init;
			return false;
		} else {
			cerr << key << " required" << endl;
			exit(1);
		}
	}
}

bool data::get_string (const char *k, string &s, bool required) {
	string key(k);
	vector<string> toks;
	functions::tokenize(k,toks,"|");
	for (int i  = 0 ; i < toks.size(); i++){
		io::println ("tok =  " + toks[i],4);
		if (configMap.find (toks[i]) != configMap.end()){
			s = configMap[toks[i]];
			return true;
		}
	}
	if (!required)
		return false;
	else {
		cerr << key << " required" << endl;
		exit(1);
	}
}
/*
int main (int argc, char *argv[]) {
	static const char *optString = "p:m:Il:o:vh?";

	static const struct option longOpts[] = {
		{ "no-index", no_argument, NULL, 'I' },
		{ "language", required_argument, NULL, 'l' },
		{ "output", required_argument, NULL, 'o' },
		{ "verbose", no_argument, NULL, 'v' },
		{ "randomize", no_argument, NULL, 0 },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, no_argument, NULL, 0 }
	};
	io::debug = 1;
	data *d = new data (argc, argv, optString, longOpts);
}
*/
