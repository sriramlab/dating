#ifndef GENOTYPE_0_H

#define GENOTYPE_0_H

#include "packedgtype.h"
#include "std.h"
#include "snp.h"
#include "printable.h"
#include "data.h"
#include "ind.h"


class genotype : public printable{
	public:
	int nsnps;
	int nchr;
	int nind;

	string genofile;
	string weightfile;
	bool givenweights;
	string snpfile;
	string indivfile;
	string outfile;
	bool outputfail;
	gsl_rng *r;

	data *d;
	
	unordered_map<string, double> maxppos;
	unordered_map<string, double> minppos;
	unordered_map<string, double> maxgpos;
	unordered_map<string, double> mingpos;
	unordered_map<string, vector<snp> > snps;       		
	unordered_map<string, int> numsnps;
    vector<unsigned int> cumnumsnps;
	unordered_set<snp,boost::hash<snp> > allsnps ;
	unordered_map<string, pair<string,int> > allsnpsmap;

    bool pack;
    packedgtype *pg;

	bool isgenotype;
	bool givenind ;
	vector<string> chrs;

	unordered_map <string, double> weightmap;
    unordered_map <string, int> indivmap;
	vector<ind> indiv;

//	genotype (string snpfile, string indivfile, string genofile, bool isgenotype=true) ; 
	genotype (string snpfile, string genofile, bool isgenotype=true, bool pack = false) ; 
    genotype (string snpfile, string indivfile, string genofile, string format = "eigenstrat", bool isgenotype = true, bool pack = false) ;
    void read_vcf (string genofile, bool isgenoytpe) ;
	genotype (data *);
	void set_weights (string weightfile); 
	void read_snps (string filename);
	void read_genotypes (string filename);
	void read_ind (string filename) ;
	void read_weights (string filename) ; 
	int get_index (string id);
	pair<string,int> get_snp (string id);
	int get_geno (string chr, int snpindex, int ind);
    int get_geno (int chrindex, int snpindex, int ind);
    inline int operator () (int chrindex, int snpindex, int ind);
	string to_string () const;

};


inline int genotype::get_geno (string chr, int snpindex, int ind){
    	snp &s = snps[chr][snpindex];
	    return s.gtype[ind];
}

inline int genotype::get_geno (int chrindex, int snpindex, int ind){
    if (!pack) { 
    	snp &s = snps[chrs[chrindex]][snpindex];
	    return s.gtype[ind];
    } else {
        pg->get_geno(chrindex, snpindex, ind);
    }
}

inline int genotype::operator () (int chrindex, int snpindex, int ind){
    if (!pack) { 
    	snp &s = snps[chrs[chrindex]][snpindex];
	    return s.gtype[ind];
    } else {
        (*pg)(chrindex, snpindex, ind);
    }
}



#endif
