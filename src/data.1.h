#ifndef DATA_1_H

#define DATA_1_H

#include "std.h"

class data {

	public:
	data (string filename);
	data (int argc, char *argv[], const char *optString, const option *longOpts) ;

	void read_config (string filename);
	void config_simple () ;
	bool get_boolean (const char *k, bool &v,  bool init, bool required = false) ;
	bool get_int (const char *k, int &v, int init, bool required = false) ;
	bool get_double (const char *k, double &v, double init, bool required = false) ;
	bool get_string (const char *k, string &s, bool required = false) ;


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
