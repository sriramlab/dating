#include "data.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"

ofstream outfs;
ofstream pairfs;


void computedpairs (int argc, char *argv[]) {
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
	
	p->get_double ("maxdis",maxdis,0.01);
	p->get_double ("binsize",binsize,0.01);
	string outfile;
	out = p->get_string ("output",outfile,true);
	if (out) {
	       outfs.open (outfile.c_str());	
	}




}

 
