#include "data.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "mergeit.h"

ofstream outfs;

mergeit::mergeit (int argc, char *argv[]) {
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
    g1 = new genotype (snp1, ind1, geno1, "eigenstrat");


	p->get_boolean ("addsnps", addsnps, false, false);
	p->get_boolean ("keepgeneticdist",keepgen,false,false);
	string geno2, snp2, ind2;
//	p->get_string ("geno2",geno2,true);
	p->get_string ("snp2",snp2,true);
//	p->get_string ("ind2",ind2,true);
//	g2 = new genotype (snp2, ind2, geno2);
	g2  = new snpmap (snp2);

}

void mergeit::merge ()  {

	p->get_string("snpoutname|snpoutfilename",snpoutfile,true);
	io::println ("snpoutfile =  " + snpoutfile, 2);
	p->get_string("indivoutname|indoutfilename",indivoutfile,true);
	p->get_string("genotypeoutname|genooutfilename",genotypeoutfile,true);


	ofstream snpfs (snpoutfile.c_str());
	ofstream indivfs (indivoutfile.c_str());
	ofstream genofs (genotypeoutfile.c_str());


	int nind = g1->nind;
	int nchr = g1->nchr;
	for (int i = 0 ; i < nchr; i++){  
		vector<snp> &tmpsnps =  g1->snps[g1->chrs[i]];
		
		vector<snp> &newsnps = g2->snps[g1->chrs[i]];
		int k = 0;
		for (int j  = 0 ; j < tmpsnps.size();j++) {
			snp &oldsnp = tmpsnps[j];
			while (k < newsnps.size() && newsnps[k].physpos<oldsnp.physpos) {
				snp &newsnp = newsnps[k];
				if (keepgen)
					snpfs << newsnp.id << "\t" << newsnp.chr <<"\t"<<newsnp.genpos<<"\t"<<newsnp.getphyspos() <<"\t"<<newsnp.var<<"\t"<<newsnp.ref<<endl;
				else {
				       double lp = j==0?0:tmpsnps[j-1].physpos;	
				       double rp = oldsnp.physpos;	
				       // Set genetic distances to the left of the leftmost SNP to zero
				       double lg = j==0?tmpsnps[j].genpos:tmpsnps[j-1].genpos;
				       double rg = oldsnp.genpos;

				       double newg = ((newsnp.physpos - lp)/(rp-lp))*(rg-lg) + lg;
					snpfs << newsnp.id << "\t" << newsnp.chr <<"\t"<<newg<<"\t"<<newsnp.getphyspos() <<"\t"<<newsnp.var<<"\t"<<newsnp.ref<<endl;
				}

				for (int l = 0 ; l < nind; l++)
					genofs << "0";
				genofs<<endl;
				k++;
			}
			// If there are duplicates, chose the current set
			if (k < newsnps.size() && newsnps[k].physpos==tmpsnps[j].physpos)
				k++;
			snpfs << oldsnp.id << "\t" << oldsnp.chr <<"\t"<<oldsnp.genpos<<"\t"<<oldsnp.getphyspos() <<"\t"<<oldsnp.var<<"\t"<<oldsnp.ref<<endl;
			vector<int> &gtype = tmpsnps[j].gtype;
			for (int l = 0 ; l < nind; l++)
				genofs << gtype[l];
			genofs<<endl;
		}
		for (; k < newsnps.size(); k++) {
			snp &newsnp = newsnps[k];
			if (keepgen)
				snpfs << newsnp.id << "\t" << newsnp.chr <<"\t"<<newsnp.genpos<<"\t"<<newsnp.physpos <<"\t"<<newsnp.var<<"\t"<<newsnp.ref<<endl;
			else { 
				snpfs << newsnp.id << "\t" << newsnp.chr <<"\t"<<tmpsnps[tmpsnps.size()-1].genpos<<"\t"<<newsnp.physpos <<"\t"<<newsnp.var<<"\t"<<newsnp.ref<<endl;
			}

			for (int l = 0 ; l < nind; l++)
				genofs << "0";
			genofs<<endl;
		}
	}
	snpfs.close();
	genofs.close();
	for (int l = 0 ; l < nind; l++){
		ind &i = g1->indiv[l];
		indivfs << i.id <<"\t" <<i.gender <<"\t"<<i.group<<endl; 
	}
	indivfs.close();
}


int main (int argc, char* argv[]) {
	mergeit m = mergeit(argc,argv);
	m.merge ();
}
