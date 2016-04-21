#include "std.h"
#include "data.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"


class physicalgenetic { 
	public:
		physicalgenetic (int argc, char *argv[]) ;
		void test ();

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

		snpmap *g;
		double maxdis;
		double binsize;
		int nbins;
};

physicalgenetic::physicalgenetic (int argc, char *argv[]) {
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
	p->get_string ("snpname",snp1,true);
	g = new snpmap (snp1);
	p->get_double ("maxdis",maxdis,1e6);
	p->get_double ("binsize",binsize,1e3);
}

void physicalgenetic::test ()  {
	nbins =  mathfn::round (maxdis/binsize );
	int nchr = g->nchr;
	vector<double>  rate(nbins);
	vector<double>  mean(nbins);
	vector<double>  sd(nbins);
	vector<int>  denom(nbins);
	double total = 0;
	double tdenom = 0 ;
	for (int i = 0 ; i < nchr; i++){  
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];

		for (int j  = 0 ; j < tmpsnps.size();j++) {
			for (int k  = j+1; k < tmpsnps.size (); k++) {
				snp &s1 = tmpsnps[j];
				snp &s2 = tmpsnps[k];

				double dis = s2.genpos - s1.genpos;
				double pdis = s2.physpos - s1.physpos;
				if (s2.physpos  - s1.physpos > maxdis)
					continue;
				int bini  = (int)(pdis/binsize);
				double r = 0;
				if (pdis > 0 )
					r = dis/pdis;
				rate[bini]+=r;
				mean[bini] += dis;
				sd[bini] += dis*dis;
				denom[bini]++;
				total += r;
				tdenom ++;

			}
		}
	}
	if (tdenom  > 0)
		total/=tdenom;
	cout << "# " << total << "\t" << tdenom << endl;
	for (int i = 0 ;i  < nbins; i++) {
		if (denom[i]>0) { 
			rate[i]/=denom[i];
			mean[i] /= denom[i];
			sd[i] /= denom[i];
			sd[i] -= mean[i]*mean[i];
		}
		cout << i << "\t" << rate[i] << "\t" << mean[i] << "\t" << sd[i] << "\t" << denom[i] << endl;
	}
}


int main (int argc, char* argv[]) {
	io::debug = 2;
	physicalgenetic c = physicalgenetic(argc,argv);
	c.test ();
}
