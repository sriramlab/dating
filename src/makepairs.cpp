#include "data.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snp.h"
#include "snpmap.h"
#include<boost/filesystem/operations.hpp>
namespace bf = boost::filesystem;		 //create an alias

ofstream outfs;
ofstream pairfs;

class makepairs { 
	public:
		makepairs (int argc, char *argv[]) ;

		data *d;
		data *p;
		snpmap *g;


		double maxdis;
		bool outputpair;
		string outputpairfile;
	
};

makepairs::makepairs (int argc, char *argv[]) {
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
	p->get_int ("debug", io::debug, 0);
	string snpfile;
	p->get_string ("snpname", snpfile, true);
	g = new snpmap (snpfile);
	
	outputpair = p->get_string ("outputpair",outputpairfile,true);
	if (outputpair) { 
		pairfs.open (outputpairfile.c_str());
	}
	int nchr = g->nchr;
	for (int i = 0 ; i < nchr; i++){  
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];
		for (int j  = 0 ; j < tmpsnps.size();j++) {
			for (int k  = j+1; k < tmpsnps.size (); k++) {
				snp &s1 = tmpsnps[j];
				snp &s2 = tmpsnps[k];
				double dis = s2.genpos - s1.genpos;
				double pdis = s2.physpos - s1.physpos;
				if (s2.genpos  - s1.genpos > maxdis)
					continue;
				if (outputpair) { 
					pairfs <<  s1.chr << "\t" <<  s1.getphyspos() << "\t" << s2.chr << "\t" << s2.getphyspos() << "\t" << dis << "\t" << pdis << endl;
				}


			}
		}
	}
	if (outputpair)
		pairfs.close ();
}





int main (int argc, char *argv[]) {
	makepairs cd (argc, argv);
}
