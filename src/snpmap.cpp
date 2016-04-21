#include "snpmap.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"



snpmap::snpmap (string snpfile) { 
//	read_test (snpfile);
	read_snps (snpfile);
}


void snpmap::read_snps (string filename){
	io::println ("Reading SNPs",2);
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	vector<snp> tmp;
	int j = 0 ;
	maxgpos = unordered_map<string, double> ();
	mingpos = unordered_map<string, double>();
	maxppos = unordered_map<string, double> ();
	minppos = unordered_map<string, double>();
	chrs = vector<string> ();
	double prevgpos ;
	double prevppos;
	string prevchr="";
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);
		
		if (line.empty())
			continue;

		string id;
		string chr;
		double genpos;
		double physpos;
		char var;
		char ref;

		ss >> id;
		ss >> chr;
		ss >> genpos;
		ss >> physpos;
		ss >> var;
		ss >> ref;

		if (chr.compare(prevchr)==0) {
			if (physpos >= prevppos){
				prevgpos = genpos;
				prevppos = physpos;
			} else {
				cerr << "Invalid snp file " << filename <<" : Line " << linenum  << "\t Must have non-decreasing physical positions" << endl;
				exit (1);
			}
		} else {
			prevchr = chr;
			prevgpos = genpos;
			prevppos = physpos;
			if (snps.find(chr) != snps.end ()) {
				cerr << "Invalid snp file "<<filename << " : Line " << linenum << "\t SNPs on a chromosome must be consecutive" << endl;
				exit (1);
			}
			chrs.push_back(chr);
			snps[chr] = vector<snp>();
		}

		snps[chr].push_back (snp(id,chr,physpos,genpos,var,ref));
		j++;

		if ( maxppos.find(chr) == maxppos.end()) {
			maxppos[chr] = physpos;
		} else if (physpos > maxppos[chr])
			maxppos[chr] = physpos;


		if ( minppos.find(chr) == minppos.end()) {
			minppos[chr] = physpos;
		} else if (physpos < minppos[chr])
			minppos[chr] = physpos;


		if ( maxgpos.find(chr) == maxgpos.end()) {
			maxgpos[chr] = genpos;
		} else if (genpos > maxgpos[chr])
			maxgpos[chr] = genpos;


		if ( mingpos.find(chr) == mingpos.end()) {
			mingpos[chr] = genpos;
		} else if (genpos < mingpos[chr])
			mingpos[chr] = genpos;
	}
	nsnps = j;
	nchr = chrs.size();

	for (int i  = 0 ; i < nchr; i++) { 
		numsnps[chrs[i]] = snps[chrs[i]].size();
		allsnps.insert ( snps[chrs[i]].begin(),snps[chrs[i]].end());
		vector<snp> &tmpsnps = snps[chrs[i]];
		for (int j = 0 ; j < tmpsnps.size(); j++){ 
			snp &s= tmpsnps[j];
			if (allsnpsmap.find(s.id) != allsnpsmap.end()){
				cerr << "Duplicate SNP ids " << s.id <<endl;
				exit(1);
			}
			allsnpsmap[s.id] = pair<string,int>(s.chr,j);
		}
	}

	inp.close();
	
	io::println ("Finished reading SNPs",2);
}


string snpmap::to_string () const {
	string s = "";
	for (int i  = 0 ; i < nchr; i++) { 
		unordered_map<string, vector<snp> >::const_iterator it  = snps.find(chrs[i]);
		const vector<snp>& tmpsnps = it->second;
		
		for (int j = 0 ; j < tmpsnps.size(); j++) {
			s = s +  tmpsnps[j].to_string() +  "\n";
		}
	}
	return s;
}


void snpmap::read_test (string filename){
	io::println ("Reading SNPs",0);
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	vector<string> tmp;
	int j = 0 ;
	int linenum  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		
		if (line.empty())
			continue;

		tmp.push_back (line);
	}
	inp.close ();
	io::println ("Finished reading SNPs",0);
}

/*
int main (int argc, char *argv[]) { 
	io::debug = 2;
	snpmap s (argv[1]);
}*/
