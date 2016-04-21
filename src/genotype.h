#ifndef GENOTYPE_H

#define GENOTYPE_H

#include "packedgtype.h"
#include "std.h"
#include "snp.h"
#include "printable.h"
#include "data.h"
#include "ind.h"
#include "intervalmap.h"


class genotype : public printable{
	public:
	int nsnp;
    int originalnsnp;
	int nchr;
	int nind;
    int originalnind;

	string genofile;
	string weightfile;
	bool givenweights;
	string snpfile;
	string indivfile;
	string outfile;
    string goodsnpfile;
    bool givengoodsnp;
    bool givenbadsnp;
    string indivfilterfile;
    bool givenindivfilter;
    int indivtofilter;
    int snptofilter;

	bool outputfail;
	gsl_rng *r;

	data *d;
	
	unordered_map<string, double> maxppos;
	unordered_map<string, double> minppos;
	unordered_map<string, double> maxgpos;
	unordered_map<string, double> mingpos;
	unordered_map<string, vector<snp> > snps;       		
    // Number of SNPs on each chromosome
	unordered_map<string, int> numsnps;
    vector<unsigned int> cumnumsnps;
	unordered_set<snp,boost::hash<snp> > allsnps ;
	unordered_map<string, pair<string,int> > allsnpsmap;
    vector<int> prevgtype;

    bool pack;
    packedgtype *pg;
    packedgtype *anc;

	bool isgenotype;
    bool hasmissing;
	bool givenind ;

    // Names of chromosomes
	vector<string> chrs;
    // Name to index map
    unordered_map <string, int> chrindexmap;

	unordered_map <string, double> weightmap;
    unordered_map <string, int> indivmap;
	vector<ind> indiv;
    unordered_map <string, vector<int > > popindexmap;
    // True => exclude
    vector<bool> *indivfilter;
    vector<bool> *snpfilter; 

    string sep;
    int prevsnpindex;

	genotype (string snpfile, string genofile, bool isgenotype=true, bool hasmissing = true, bool pack = true, string sep = "" ) ; 
    genotype (string snpfile, string indivfile, string genofile, string format = "eigenstrat", bool isgenotype = true, string indivfilterfile = "", string goodsnpname = "", bool hasmissing = true, bool pack = true, string sep = "") ;
    void read_vcf (string genofile, bool isgenoytpe, intervalmap *imap, string goodsnpsname, bool isx = false) ;
    void read_vcf (string genofile, bool isgenotype = true, string filterpolicy = "", string goodsnpsname =  "", bool isx = false) ;
	genotype (data *);
	void set_weights (string weightfile); 
	void read_snps (string filename, string goodsnpname = "");
    void read_packedes (string filename) ;
	void read_es (string filename);
	void read_genotypes (string filename);
	void read_ind (string filename, string filterfile) ;
	void read_weights (string filename) ; 
	int get_index (string id);
	pair<string,int> get_snp (string id);

	int get_geno (string chr, unsigned int snpindex, int ind);
    int get_geno (int chrindex, unsigned int snpindex, int ind);
    int get_geno (unsigned int snpindex, int ind);
    int get_geno (string chr, unsigned int snpindex, vector <int> &v) ;
    int get_geno (int chrindex, unsigned int snpindex, vector<int> &v) ;
    int get_geno (unsigned int snpindex, vector<int> &v) ;



    int get_prob (string chr, unsigned int snpindex, vector<double> &v) ;
    int get_prob (int chrindex, unsigned int snpindex, vector<double> &v) ;
    int get_prob (unsigned int snpindex, vector<double> &v);

    void set_geno (int chrindex, unsigned int snpindex, unsigned int ind, int val) ;
    void set_geno (int chrindex, unsigned int snpindex, vector<int> &val) ;
    void write_genotypes(string filename) ;
    void set_freq ( );
    double get_freq (vector<int> &gtype, string popname = "") ;
    double get_freq (vector<int> &gtype, int &vc, int &rc, int &mc, string popname = "") ;
    double get_freq (int chrindex, unsigned int snpindex, int &vc, int &rc, int &mc, string popname = "") ;
    double get_freq (unsigned int snpindex, int &vc, int &rc, int &mc, string popname = "") ; 
    double get_freq (unsigned int snpindex, string popname = "") ;
    double get_freq (int chrindex, unsigned int snpindex, string popname = "") ;

    ostream& write_genotypes ( ostream & os) ;

    inline int operator () (int chrindex, int snpindex, int ind);
	string to_string () const;

};


inline int genotype::get_geno (string chr, unsigned int snpindex, int ind){
    if (!pack) { 
    	snp &s = snps[chr][snpindex];
	    return s.gtype[ind];
    } else { 
        pg->get_geno(chrindexmap[chr], snpindex, ind);
    }
}

inline int genotype::get_geno (int chrindex, unsigned int snpindex, int ind){
    if (!pack) { 
    	snp &s = snps[chrs[chrindex]][snpindex];
	    return s.gtype[ind];
    } else {
        pg->get_geno(chrindex, snpindex, ind);
    }
}

inline int genotype::get_geno (unsigned int snpindex, int ind){
    if (!pack) { 
        int i = 0 ;
        for (; i < cumnumsnps.size() && cumnumsnps[i] <= snpindex; i++);
    	snp &s = snps[chrs[i-1]][snpindex];
	    return s.gtype[ind];
    } else {
        pg->get_geno(snpindex, ind);
    }
}


inline int genotype::get_geno (string chr, unsigned int snpindex, vector <int> &v) {
    get_geno (chrindexmap[chr], snpindex, v);
}

inline int genotype::get_geno (int chrindex, unsigned int snpindex, vector<int> &v) { 
    if (!pack) {
    	snp &s = snps[chrs[chrindex]][snpindex];
        for (int i = 0 ; i < nind; i++)
            v[i] = s.gtype[i];
    } else { 
        pg->get_geno (chrindex, snpindex,v);
    }

}

inline int genotype::get_geno (unsigned int snpindex, vector<int> &v) { 
    if (!pack) {
        int i = 0 ;
        for (; i < cumnumsnps.size() && cumnumsnps[i] <= snpindex; i++);
    	snp &s = snps[chrs[i]][snpindex];
        for (int i = 0 ; i < nind; i++)
            v[i] = s.gtype[i];
    } else { 
        pg->get_geno ( snpindex,v);
    }

}

inline int genotype::get_prob (string chr, unsigned int snpindex, vector<double> &v) {
    get_prob (chrindexmap[chr], snpindex, v);
}

inline int genotype::get_prob (int chrindex, unsigned int snpindex, vector<double> &v) { 
    if (!pack) {
        cerr << "Error in genotype::get_prob: prob entry needs packed format" << endl;   
        exit (1);     
    } else  {
        pg->get_prob (chrindex, snpindex, v);
    }
}

inline int genotype::get_prob (unsigned int snpindex, vector<double> &v) { 
    if (!pack) {
        cerr << "Error in genotype::get_prob: prob entry needs packed format" << endl;   
        exit (1);     
    } else  {
        pg->get_prob (snpindex, v);
    }
}

inline void genotype::set_geno (int chrindex, unsigned int snpindex, unsigned int ind, int val) {
    if (!pack) {
    	snp &s = snps[chrs[chrindex]][snpindex];
        int oldval = s.gtype[ind];
        s.gtype[ind] = val;
        if (val != oldval) {
            int a = (val!=9)?val:0;
            int b = (oldval!=9)?oldval:0;
            int vc = s.vcount + a - b;
            a = (val!=9)?(isgenotype?(2-val):(1-val)):0;
            b = (oldval!=9)?(isgenotype?(2-oldval):(1-oldval)):0;
            int rc = s.rcount + b - a;

            s.vcount = vc;
            s.rcount = rc;
            s.freq = (vc + rc > 0)?((1.*vc)/(vc+rc)):0;
        }

    } else {
        pg->set_geno (chrindex, snpindex, ind, val);
    }
}

inline void genotype::set_geno (int chrindex, unsigned int snpindex, vector<int> &val) {
    if (!pack) {
        int vc = 0;
        int rc = 0;
        snp &s = snps[chrs[chrindex]][snpindex];
        for (int i = 0 ; i < nind ; i++) {
            int oldval = s.gtype[i];
            s.gtype[i] = val[i];
            if (val[i]!=9) {
                vc += val[i];
                rc += isgenotype?(2-val[i]):(1-val[i]);
            }
        }
        s.vcount = vc;
        s.rcount = rc;
        s.freq = (vc + rc > 0)?((1.*vc)/(vc+rc)):0;
    } else {
        pg->set_geno (chrindex, snpindex, val);
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
