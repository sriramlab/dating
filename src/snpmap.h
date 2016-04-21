#ifndef SNPMAP_H

#define SNPMAP_H

#include "std.h"
#include "snp.h"
#include "printable.h"
#include "data.h"
#include "ind.h"


class snpmap : public printable{
	public:
	int nsnps;
	int nchr;

	
	unordered_map<string, double> maxppos;
	unordered_map<string, double> minppos;
	unordered_map<string, double> maxgpos;
	unordered_map<string, double> mingpos;
	unordered_map<string, vector<snp> > snps;       		
	unordered_map<string, int> numsnps;
	unordered_set<snp,boost::hash<snp> > allsnps ;
	unordered_map<string, pair<string,int> > allsnpsmap;

	vector<string> chrs;

	snpmap (string snpfile) ; 
	void read_snps (string filename);
	void read_test (string filename);
	string to_string () const;

};
#endif
