#include "genotype.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"


genotype::genotype (data *d) { 
	this->d = d;
	d->get_string ("genotypename",genofile,true);
	d->get_string ("snpname",snpfile,true);
	d->get_string ("indivname",indivfile,true);
	d->get_boolean ("genotype", isgenotype, true);

	read_snps (snpfile);
	read_ind (indivfile);
	read_genotypes (genofile);

	givenweights = d->get_string ("weightname",weightfile,false);
	if (givenweights) 
		read_weights (weightfile);
}

genotype::genotype (string snpfile, string indivfile, string genofile, bool isgenotype) { 
	read_snps (snpfile);
	read_ind (indivfile);
	this->isgenotype = isgenotype;
	read_genotypes (genofile);
}

genotype::genotype (string snpfile, string genofile, bool isgenotype) { 
	read_snps (snpfile);
	this->isgenotype = isgenotype;
	givenind = false;
	read_genotypes (genofile);
}

void genotype::set_weights (string weightfile){ 
	givenweights  = true;
	read_weights (weightfile);
}

void genotype::read_genotypes (string filename) { 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	int chrindex = 0;
	int snpindex = 0;
	bool setind = false;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		
		if (line.empty())
			continue;

		vector<snp> &sv = snps[chrs[chrindex]];
		vector<int> &gtype = sv[snpindex].gtype;
		io::println ("Length = "  + tostring(line.length()),2);
		int vc = 0;
		int rc = 0;

		if (!givenind ) { 
			if (!setind){
				nind = line.length();
				for (int i  =0 ;i  < nind; i++){
					string id = tostring(i);
					indiv.push_back (ind (id, 'U', "G"));
				}
				setind = true;
			} else {
				if (nind != line.length()) {
					cerr << "Bad genotype in " << filename << ": Line "<< linenum  << endl;
					exit(1);
				}
			}
			gtype.resize(nind,9);
		} else { 
			if (nind != line.length()) {
				cerr << "Bad genotype in " << filename << ": Line "<< linenum  << endl;
				exit(1);
			}
		}

		for (int i = 0 ; i < line.length();i++){
			

			char c = line.at(i);
//			gtype[i] = atoi (&c);
			gtype[i] = c-48;
			if (gtype[i]==9 || (gtype[i]>=0 && gtype[i]<=2)){
				if (!isgenotype && gtype[i]==2) {
					cerr << "Data is not haploid in " << filename << ": Line "<< linenum <<", col = " << i << endl;
				}
				if (gtype[i]!=9) {
					vc += gtype[i];
					rc = isgenotype?2:1;
					rc = rc - gtype[i];
				}
			} else {
				cerr << "Bad genotype in " << filename << ": Line "<< linenum <<", col = " << i << endl;
				exit(1);
			}
			if ( io::debug >= 3) {
				cout << "gtype " << tostring(i) << "  = " << tostring (c) << "," << tostring (gtype[i])<<endl;
			}
		}
		sv[snpindex].vcount = vc;
		sv[snpindex].rcount = rc;
		sv[snpindex].freq = vc+rc==0? 0:((double)vc)/(vc+rc);

		snpindex ++;	
		if (snpindex < sv.size()){}
		else {
			snpindex = 0 ;
			chrindex++;
		}
		
	}	
	inp.close();

	if (io::debug >= 2) { 
		for (int i = 0  ; i < nchr; i++){
			vector<snp> &tmpsnps = snps[chrs[i]];
			for (int j =  0 ; j < tmpsnps.size();j++){
				cout << tmpsnps[j].getgtype () << endl;
			}
			cout << endl;
		}	
	}
}

void genotype::read_weights (string filename) { 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
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
		double w;

		ss >>id;
		ss >>w;
		weightmap[id] = w;
	}
	inp.close ();

}

void genotype::read_ind (string filename) { 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	vector<snp> tmp;
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
		char g;
		string group;

		ss >>id;
		ss >>g;
		ss >>group;
		ind x (id, g, group);
		indiv.push_back (ind (id,g,group));
	}
	nind  = linenum;
	inp.close ();

	for (int i = 0 ; i < nchr; i++){
		vector<snp> &tmpsnps =  snps[chrs[i]];
		for (int  j = 0 ; j < tmpsnps.size(); j++){
			tmpsnps[j].gtype.resize(nind,9);
		}
	}
	givenind = true;
}



void genotype::read_snps (string filename){
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
	int totallines = 0;
	int snpcount = 0 ;
	
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
		string var;
		string ref;

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
			if (prevchr.compare ("") !=0){
				snps[prevchr] = vector<snp>(snpcount);
			}
			prevchr = chr;
			prevgpos = genpos;
			prevppos = physpos;
			if (snps.find(chr) != snps.end ()) {
				cerr << "Invalid snp file "<<filename << " : Line " << linenum << "\t SNPs on a chromosome must be consecutive" << endl;
				exit (1);
			}
			chrs.push_back(chr);
			snpcount = 0;
		}

		totallines++;
		snpcount++;

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
	snps[prevchr] = vector<snp>(snpcount);

	/*
	for (int i = 0 ; i < chrs.size(); i++){
		cout << i << "\t" <<  chrs[i] << "\t" << snps[chrs[i]].size() << endl;
	}
	*/

	inp.close ();
	inp.open (filename.c_str());


	linenum = 0;
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
		string var;
		string ref;

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
			snpcount ++;
		} else {
			prevchr = chr;
			prevgpos = genpos;
			prevppos = physpos;
			snpcount =  0;
		}

		snps[chr][snpcount] = snp(id,chr,physpos,genpos,var,ref);
		j++;

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
	
}

string genotype::to_string () const {
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

int genotype::get_index (string id) { 
	if (allsnpsmap.find (id) == allsnpsmap.end())
		return -1;
	else {
		pair<string, int> p = allsnpsmap[id];
		return p.second;
	}
}

pair<string,int> genotype::get_snp (string id) { 
	pair<string,int>  p ("",-1);
	if (allsnpsmap.find (id) != allsnpsmap.end())
		p = allsnpsmap[id];
	return p;
}
