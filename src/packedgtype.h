#ifndef PACKEDGTYPE_H

#define PACKEDGTYPE_H

#include "std.h"
#include "snp.h"
#include "printable.h"
#include "data.h"
#include "ind.h"


class packedgtype  {
    public:
    unsigned char *gtype;
    unsigned char *tgtype;
    int *vcount;
    int *rcount;
    unsigned char mask;
    int wordsize;
    unsigned int unitsperword;
    int unitsize;
    int order;
    int nind;
    int nsnp;
    unsigned int ncol;
    unsigned int nrow;
    bool isgenotype;
    bool hasmissing;
    bool isdouble;
    string entryformat ;
//    typedef enum {GENO, PROB} etype;
    typedef int etype;
    etype entrytype;

    vector<unsigned int> cumnumsnps;
    vector<snp> *snps;
    data *d;

    static const int SNPMAJOR = 0;
    static const int INDMAJOR = 1;

    static const int GENO = 0;
    static const int PROB = 1;

    packedgtype (int nind, int nsnp, vector<unsigned int> cumnumsnps, bool isgenotype = true, bool hasmissing = true, data *d = NULL, vector<snp> *s = NULL, string entryformat = "geno", int order = SNPMAJOR) ;
    packedgtype (int nsnp, vector<unsigned int> cumnumsnps, bool isgenotype = true, bool hasmissing = true , data *d = NULL, vector<snp> *s = NULL, string entryformat = "geno", int order = SNPMAJOR) ;

    int get_geno (unsigned int snp, unsigned int ind) ;
    int get_geno (int chr, unsigned int snp, unsigned int ind) ;
    void get_geno (int chr, unsigned int snp, vector<int> &v) ; 
    void get_geno (unsigned int snp, vector<int> &v);


    double get_prob (unsigned int snp, unsigned int ind) ;
    double get_prob (int chr, unsigned int snp, unsigned int ind) ; 
    void get_prob (int chr, unsigned int snp, vector<double> &v) ; 
    void get_prob (unsigned int snp, vector<double> &v);

    void set_geno (int chr, unsigned int snp, unsigned int ind, int val) ;
    void set_geno (int chr, unsigned int snp, vector<int> &val) ;
    void set_geno (unsigned int snp, vector<int> &val) ;
    
    int  operator () (int chr, unsigned int snp, unsigned int ind) ;
    int  operator () (unsigned int snp, unsigned int ind) ;

    void set_freq ( int chr, vector<snp> &s) ;
    int get_vc (int chr, unsigned int snp) ;
    double get (int chr, unsigned int snp, unsigned int ind) ;
    double get ( unsigned int snp, unsigned int ind) ;
    
    etype format_to_type (string &) ;
    string get_format (int e) ; 
    void set_metadata () ;

    void write_packedes(string filename) ;
    ostream& write_packedes ( ostream &ofs) ;
    void read_packedes(string filename) ;
    void read_packedes(string filename, vector<bool> *indivfilter, vector<bool> *snpfilter) ;

    void read_es (string filename, vector<bool> *indivfilter, vector<bool> *snpfilter) ;
    void read_eigenstrat (string filename, vector<bool> *indivfilter, vector<bool> *snpfilter) ;
    static void print_packedes(string filename) ;

};

inline int packedgtype::get_geno (unsigned int snp, unsigned int ind) { 
    unsigned char c = gtype [ snp *ncol + ind/unitsperword];
    int y = (c>>unitsize*(ind%unitsperword))&mask;
    y=(y==3)?9:y;
    return (y);
}

inline int packedgtype::get_geno (int chr, unsigned int snp, unsigned int ind) { 
    return get_geno (cumnumsnps[chr] + snp, ind);
}


inline void packedgtype::set_geno (int chr, unsigned int snp, unsigned int ind, int val) {
    int  i = (cumnumsnps[chr]+snp) *ncol + ind/unitsperword;
    unsigned char c = gtype [ (cumnumsnps[chr]+snp) *ncol + ind/unitsperword];
    unsigned char d = (c>>unitsize*(ind%unitsperword))&mask;
    int oldval = d; oldval=(oldval==3)?9:oldval;

    unsigned char e = val; 
    d = d ^ val;
    d = d<<unitsize*(ind%unitsperword); 
    c=c^d; // zeros out the right character and replaces with val
    gtype [ (cumnumsnps[chr]+snp) *ncol + ind/unitsperword] = c; 
    if (oldval!=val) {
        int a = (val!=9)?val:0;
        int b = (oldval!=9)?oldval:0;
        int vc = vcount [cumnumsnps[chr]+snp ] + a - b;
        vcount[cumnumsnps[chr]+snp] = vc;
        a = (val!=9)?(isgenotype?(2-val):(1-val)):0;
        b = (oldval!=9)?(isgenotype?(2-oldval):(1-oldval)):0;
        int rc = rcount [cumnumsnps[chr]+snp ] + b - a;
        (*snps)[snp].vcount = vc;
        (*snps)[snp].rcount = rc;
        (*snps)[snp].freq = (vc + rc > 0)?((1.*vc)/(vc+rc)):0;
    }
}

inline void packedgtype::set_geno (int chr, unsigned int snp, vector<int> &val) {
    set_geno (cumnumsnps[chr]+snp, val);
}

inline void packedgtype::set_geno (unsigned int snp, vector<int> &val) {
    int  i = snp*ncol ;
    int vc = 0;
    int rc = 0;
    int shift = 0;
    for (int j = 0  ; j < nind ; j++) { 
        if (j%unitsperword==0)
            gtype[i+j/unitsperword]=0;
        unsigned char c = (char)val[j];
        c=(c==9)?3:c;
        gtype[i+j/unitsperword] = gtype[i+j/unitsperword] |(c<< (unitsize*shift)) ;
        if (c!=9) {
            vc += c;
            rc += isgenotype?(2-c):(1-c);
        }
        shift = (shift +1 ) % unitsperword;
    }
    vcount[snp] = vc;
    rcount[snp] = rc;

}


inline double packedgtype::get (int chr, unsigned int snp, unsigned int ind) { 
    return get (cumnumsnps[chr]+snp, ind);
}

inline double packedgtype::get (unsigned int snp, unsigned int ind) { 
    if (entryformat.compare ("prob")==0) {
        return get_prob (snp, ind);
    } else {
        cerr << "Undefined entry type in packedgtype::get" << endl;
        exit (1);
    }
}

inline double packedgtype::get_prob (int chr, unsigned int snp, unsigned int ind) { 
    return get_prob (cumnumsnps[chr]+snp, ind);
}

inline double packedgtype::get_prob (unsigned int snp, unsigned int ind) { 
    double y;
    unsigned char c1 = gtype[snp*ncol +  2*ind];
    unsigned char c2 = gtype[snp*ncol +  2*ind + 1];
    int n = 0 ;
    n += ((c1&(mask<<unitsize))>>unitsize); 
    n = 10*n + (c1&mask);
    n = 10*n +  ((c2&(mask<<unitsize))>>unitsize); 
    n = 10*n + (c2&mask);
    y = n/10000.0;
    int x = ((c1&(mask<<unitsize))>>unitsize);
    return y;
}

inline  void packedgtype::get_prob (unsigned int snp, vector<double> &v) { 
    for (int i = 0 ;i < nind ; i++)
        v[i]  = get_prob (snp, i);

}


inline  void packedgtype::get_prob (int chr, unsigned int snp, vector<double> &v) { 
    get_prob (cumnumsnps[chr]+snp, v);
}

inline int  packedgtype::operator () (int chr, unsigned int snp, unsigned int ind) { 
    unsigned char c = gtype [ (cumnumsnps[chr]+snp) *ncol + ind/unitsperword];
    int y = (c>>unitsize*(ind%unitsperword))&mask;
    y=(y==3)?9:y;
    return (y);
}

inline int  packedgtype::operator () (unsigned int snp, unsigned  int ind) { 
    unsigned char c = gtype [ snp *ncol + ind/unitsperword];
    int y = (c>>unitsize*(ind%unitsperword))&mask;
    y=(y==3)?9:y;
    return (y);
}



inline int packedgtype::get_vc (int chr, unsigned int snp) {
   int vc = vcount [ cumnumsnps[chr]+snp];
   return (vc);
}

inline  void packedgtype::get_geno (int chr, unsigned int snp, vector<int> &v) { 
    get_geno ( cumnumsnps[chr]+snp, v);
}

inline  void packedgtype::get_geno (unsigned int snp, vector<int> &v) { 
    for (int i = 0 ;i < nind ; i++) {
        unsigned char c = gtype [ snp *ncol + i/unitsperword];
        int y = (c>>unitsize*(i%unitsperword))&mask;
        y=(y==3)?9:y;
        v[i] = y;
    }
}


#endif

