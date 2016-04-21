#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "genotype.h"

ofstream outfs;

class conditionalafs { 
	public:
		conditionalafs (int argc, char *argv[]) ;
		void convert ();

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

		// genotype data
		genotype *g;

		bool out;
		string snpoutfile;
		string indivoutfile;
		string genotypeoutfile;

		bool uniqpos;
        bool isgenotype;
        string outputformat;
};

conditionalafs::conditionalafs (int argc, char *argv[]) {
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
    vector<string> tok;
    functions::tokenize(geno1.c_str(), tok,".");
    if (tok[tok.size()-1].compare("gz")==0){
        if (tok[tok.size()-2].compare("vcf")==0)
            inputformat = "vcf";
    }

    outputformat = "eigenstrat";
    p->get_string ("outputformat", outputformat, false);
	p->get_boolean ("isgenotype", isgenotype, true);

	g = new genotype (snp1,ind1, geno1, inputformat, isgenotype);
}

void conditionalafs::convert ()  {

	string selectfile;
    string afsoutfile;
    double window;

	unordered_map<string, string > *smapsnps = NULL;
    bool selectsnps = p->get_string ("selectsnpname",selectfile, false);
    p->get_double ("window", window, 0.00025, false);
	p->get_string("afsout",afsoutfile,true);

	if (selectsnps){
		smapsnps = new unordered_map <string, string > ();
		fileio::read_map (selectfile, smapsnps,0);
		if (io::debug >= 2)  {
			typedef unordered_map<string, string>::iterator iter;
			for ( iter i  = smapsnps->begin();  i != smapsnps->end(); i++){
				cout << i->first << "\t" << i->second << endl;
			}
		}
	}
	
	ofstream afsout (afsoutfile.c_str());
	int nind = g->nind;
	int nchr = g->nchr;

    vector<int> index;
	for (int i = 0 ; i < nchr; i++){
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];
		double prevgpos ;
		bool flag = false;
		for (int j  = 0 ; j < tmpsnps.size();j++) {
			snp &s = tmpsnps[j];
            if (selectsnps) {
                if (smapsnps->find(s.id) == smapsnps->end()) 
                    continue;
            }
            vector<int> &gtype = s.gtype;
            index.clear ();
            for (int k = 0 ; k < nind; k++)
                if (gtype[k])
                    index.push_back (k);
            int denom = index.size()-1;
            vector<int> afs (denom);
            int nss = 0 ;

            for (int k=j-1; k>=0 ; k-- ) {
                snp &t = tmpsnps[k];
                if (s.genpos - t.genpos > window)
                    break;
                int ac = 0 ;
                for (int l = 0 ; l < index.size(); l++) {
                    ac += t.gtype[index[l]];
                }
                if (ac>0 && ac <= denom) {
                    afs[ac-1] ++;
                    nss ++;
                }
            }
            for (int k = j+1; k < tmpsnps.size(); k++ ){ 
                snp &t = tmpsnps[k];
                if (t.genpos - s.genpos > window)
                    break;
                int ac = 0 ;
                for (int l = 0 ; l < index.size(); l++) {
                    ac += t.gtype[index[l]];
                }
                if (ac>0 && ac <= denom) {
                    afs[ac-1] ++;
                    nss ++;
                }
            }

            for (int i = 0 ; i < afs.size(); i++ ) {
                double f = ((double)afs[i])/nss;
                int j = i+1;
                afsout << j  << "\t" <<  afs[i] << "\t" << nss << endl;
            }
            break;
		}
	}
    afsout.close();
}


int main (int argc, char* argv[]) {
	conditionalafs c = conditionalafs(argc,argv);
	c.convert ();
}
