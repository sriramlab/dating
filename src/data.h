#ifndef DATA_1_H

#define DATA_1_H

#include "std.h"

class data {

	public:
	data (string filename);
    data (string filename, map<string, string> &cmap);
	data (int argc, char *argv[], const char *optString = NULL, const option *longOpts = NULL) ;

//	void read_config (string filename);
	void read_config (int argc, char *argv[], int optind);
	void config_simple () ;
	bool get_boolean (const char *k, bool &v,  bool init, bool required = false) ;
	bool get_int (const char *k, int &v, int init, bool required = false) ;
	bool get_double (const char *k, double &v, double init, bool required = false) ;
	bool get_double_vector (const char *k, vector<double> &v, bool required) ;
	bool get_string (const char *k, string &s, bool required = false) ;
	void macro_sub ()  ;
	void print_parameters ()  ;


	map<string, string> configMap;
	map<string, string> clMap;
	string snpfile;
	string genofile;
	string indfile;
	string outputname;

	double binsize;
	double maxdis;

	int seed;

    
};
#endif
