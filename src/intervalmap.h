#ifndef INTERVALMAP_H

#define INTERVALMAP_H

#include "std.h"
#include "printable.h"
#include "snp.h"
#include "interval.h"

class snpmap;

class intervalmap:public printable  {
		
	public:
        bool givenchrom;
        string chrom;

		// chr-> vector of intervals
		unordered_map <string, vector<interval> > intervals;
        unordered_map <string, int> intsperchr;
        
		// chr-> vector of snps
		unordered_map<string, vector<snp> > snps;       		
		// contains map and non-map snps
		// All SNPs in an interval
		unordered_map<interval, vector<int>, boost::hash<interval> > snpsperint;

		// All counts of crossovers between adjacent SNPs in an interval
		unordered_map<interval, vector<int>, boost::hash<interval> > countsperint;

		// All SNPs
		unordered_set<snp,boost::hash<snp> > allsnps;
		snpmap *smap;
		vector<string> chrs;
        unordered_set <string> chrmap;
		gsl_rng *r;

        intervalmap () {}
        intervalmap (string filename, string chrom = "") ;
		intervalmap (string filename, gsl_rng *r, snpmap *);
        void getstats (istream &inp ) ;
        void read_stream (istream &inp, gsl_rng *r, snpmap *smap ) ;
        void merge ( ) ;
        void add (const intervalmap &im ) ;
        void intersect (const intervalmap &im ) ;
		string to_string() const;
		string stats () ;
};

#endif
