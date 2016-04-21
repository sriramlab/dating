#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"

ofstream outfs;

computehet::computehet (int argc, char *argv[]) {
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
	
	string geno1, snp1, ind1;
    string inputformat = "eigenstrat";
	p->get_string ("genotypename",geno1,true);
	p->get_string ("snpname",snp1,true);
	p->get_string ("indivname",ind1,true);
	p->get_boolean ("isgenotype", isgenotype, true);
	p->get_boolean ("isx", isx, false);


}
 
