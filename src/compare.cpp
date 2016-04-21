#include "std.h"
#include "printable.h"
#include "genotype.h"
#include "data.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"

class compare { 
	public:
		compare (int argc, char *argv[]) ;
		void comp ()  ;

		// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

		// genotype data
		genotype *g1;
		genotype *g2;

		bool out;
		string snpoutfile;
		string indivoutfile;
		string genotypeoutfile;

		bool flipsnps;
		bool badsnps;
		bool addsnps;

		bool keepgen;
};

ofstream outfs;

compare::compare (int argc, char *argv[]) {
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
	p->get_string ("geno1",geno1,true);
	p->get_string ("snp1",snp1,true);
	p->get_string ("ind1",ind1,true);
	g1 = new genotype (snp1, ind1, geno1);

	string geno2, snp2, ind2;
	p->get_string ("geno2",geno2,true);
	p->get_string ("snp2",snp2,true);
	p->get_string ("ind2",ind2,true);
	g2 = new genotype (snp2, ind2, geno2);

}

void compare::comp ()  {
	int nind = g1->nind;
	int nchr = g1->nchr;
	char c1[2*nind];
	char c2[2*nind];
	char a1[2];
	char a2[2];
	int b[3][3];
	int match = 0 ;
	int denom = 0 ;
	for (int i =0; i < 3;i++)
		for (int j = 0 ; j < 3; j++)
			b[i][j] = 0;
	for (int i = 0 ; i < nchr; i++){  
		int match1 = 0 ;
		int denom1 = 0 ;
		vector<snp> &s1 =  g1->snps[g1->chrs[i]];

		vector<snp> &s2 = g2->snps[g1->chrs[i]];
		int k = 0;
		for (int j  = 0 ; j < s1.size();j++) {
				int match2 = 0 ;
				int denom2 = 0 ;
				snp &t1  = s1[j];
				snp &t2 = s2[j];
				vector<int> &gtype1 = t1.gtype;
				vector<int> &gtype2 = t2.gtype;

				a1[1] = t1.var;
				a1[0] = t1.ref;

				a2[1] = t2.var;
				a2[0] = t2.ref;
				for (int l =  0 ; l < nind; l++) {
					if (gtype1[l]==9)
							c1[2*l] = c1[2*l+1] = 0;
					else {
						switch (gtype1[l]) {
							case 0: c1[2*l]  = c1[2*l+1] = a1[0];break;
							case 1: c1[2*l]  = a1[0]; c1[2*l+1] = a1[1];break;
							case 2: c1[2*l] = c1[2*l+1] = a1[1];break;
						}

						if (c1[2*l]>c1[2*l+1]){
							char c = c1[2*l];
							c1[2*l] =  c1[2*l+1];
							c1[2*l+1] =  c;
						}
					}
				}

				for (int l =  0 ; l < nind; l++) {
					if (gtype2[l]==9)
							c2[2*l] = c2[2*l+1] = 0;
					else  {
						switch (gtype2[l]) {
							case 0: c2[2*l]  = c2[2*l+1] = a2[0];break;
							case 1: c2[2*l]  = a2[0]; c2[2*l+1] = a2[1];break;
							case 2: c2[2*l] = c2[2*l+1] = a2[1];break;
						}

						if (c2[2*l]>c2[2*l+1]){
							char c = c2[2*l];
							c2[2*l] =  c2[2*l+1];
							c2[2*l+1] =  c;
						}
					}
				}


				for (int l =  0 ; l < nind; l++){
						if (c1[2*l]!=0 && c2[2*l] !=0 && c1[2*l+1]!=0 && c2[2*l+1]!=0) { 
							denom1 ++;
							match1 += (c1[2*l]==c2[2*l]) && (c1[2*l+1]==c2[2*l+1]);
							denom2 ++;
							match2 += (c1[2*l]==c2[2*l]) && (c1[2*l+1]==c2[2*l+1]);

							int i1 = (c1[2*l]==a2[1])+(c1[2*l+1]==a2[1]);
							int i2 = (c2[2*l]==a2[1])+(c2[2*l+1]==a2[1]);
							b[i1][i2]++;
						}
				}
				double x =  ((double)match2)/denom2;
				cout << t1.id << "\t" << t1.chr << "\t" << t1.physpos << "\t" << x << "\t" << denom2  << "\t" << match2 << endl;	
		}
		match += match1;
		denom += denom1;
		double x =  ((double)match1)/denom1;

		cout << g1->chrs[i] << "\t" << x << endl;
	}
	double x =  ((double)match)/denom;
	cout << x << endl;
	int n1 = 0;
	int n2 = 0;
	for (int i  = 0 ; i < 3; i++){
		for (int j =  0; j < 3; j++) {
				cout << b[i][j] <<"\t";
				if (i==j)
						n1+=b[i][j];
				n2+=b[i][j];
		}
		cout << endl;
	}
	x = ((double )n1)/n2;
	cout << x << endl;
	for (int j =  0; j < 3; j++) {
		int n2 = 0 ;
		for (int i  = 0 ; i < 3; i++){
				n2 += b[i][j];
		}
		x = ((double)b[j][j])/n2;
		cout << x << "\t";
	}
	cout << endl;

}


int main (int argc, char* argv[]) {
	compare m = compare(argc,argv);
	m.comp ();
}
