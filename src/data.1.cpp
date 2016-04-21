#include "std.h"
#include "data.1.h"
#include "io.h"

data::data (string filename){
	read_config(filename);
}

void data::read_config (string filename){	
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
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
		io::println ("key:" + key + "=" + value, 3);
	}
	inp.close();

	map<string,string>::iterator i;
	for ( i = clMap.begin (); i!=clMap.end (); i++) {
		string k = i->first;
		string v = i->second;
		configMap[k] = v;
	}

	io::println ("***********Parameters*************",1);
	for (i = configMap.begin(); i!=configMap.end(); i++){
		string s = i->first +"\t" + i->second ;
		io::println (s,1);
	}
	io::println ("********",1);

}


data::data (int argc, char *argv[], const char *optString, const option *longOpts) { 
	int opt = 0;
	int longIndex ;
	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	while( opt != -1 ) {
		if (io::debug >= 1)
			cout << "opt = " << opt << "\t optind = " << optind << "\t" << optarg << endl;
		switch( opt ) { 
			case 'v':
				io::debug  = 2;
			break;

			
			case 'h':  
			case '?':
				cout << "Usage:"<<endl;
				break;
				
			default:
				clMap[string(longOpts[longIndex].name)] = string(optarg);
				break;
		}

		opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	}

	if  ( argc - optind > 0) { 
		read_config (string(argv[optind]));
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
	if (configMap.find (key) != configMap.end()){
		s = configMap[key];
		return true;
	} else  {
		if (!required)
			return false;
		else {
			cerr << key << " required" << endl;
			exit(1);
		}
	}

}

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
