#include "packedgtype.h"
#include "genotype.0.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"
#include "gzstream.h"
#include "std.h"
#include "functions.h"


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
	d->get_boolean ("outputfail", outputfail, false, false);
}

/*
genotype::genotype (string snpfile, string indivfile, string genofile, bool isgenotype) { 
	read_snps (snpfile);
	read_ind (indivfile);
	this->isgenotype = isgenotype;
	read_genotypes (genofile);
}*/

genotype::genotype (string snpfile, string genofile, bool isgenotype, bool pack) { 
    this->pack = pack;
    if (pack) { 
        cerr << "Need .ind file for packed format" << endl;
        exit(1);
    }
	read_snps (snpfile);
	this->isgenotype = isgenotype;
	givenind = false;
	read_genotypes (genofile);
}

genotype::genotype (string snpfile, string indivfile, string genofile, string format, bool isgenotype, bool pack) { 
    this->pack = pack;
    if (format.compare ("eigenstrat")==0) {
    	read_snps (snpfile);
        if (functions::fileExists(indivfile))
            read_ind (indivfile);
        else
            givenind = false;
	    this->isgenotype = isgenotype;
         read_genotypes (genofile);
    } else if (format.compare ("vcf")==0) { 
        if (functions::fileExists(indivfile))
            read_ind (indivfile);
        read_vcf (genofile, isgenotype);
        this->isgenotype = isgenotype;
    } else if (format.compare ("packedeigenstrat")==0) { 
        read_snps (snpfile);
        if (functions::fileExists(indivfile))
            read_ind (indivfile);
        else
            givenind = false;
	    this->isgenotype = isgenotype;

        if (pack && !isgenotype) {
            pg = new packedgtype (nind, nsnps, cumnumsnps, isgenotype);
            pg->read_eigenstrat (genofile, NULL, NULL);
        } else  {
            cerr << "Bad format " << endl;
            exit(1);
        }
    }
}


void genotype::read_vcf (string genofile, bool isgenoytpe) { 
    igzstream ig (genofile.c_str());
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
    int size =  10000;
    int inc = 10000;
	
    int linenum = 0 ;
    bool setind = false;
    int numhap = 2;
    if (isgenotype)
        numhap = 1;
    while (std::getline (ig, line) ) {
        linenum ++;
        char c = line[0];
        if (c=='#') { 
            if (line.find ("CHROM")!=string::npos) {
                if (!givenind)  {
                    vector<string> toks;
                    functions::tokenize (line.c_str(), toks, " \t");
                    nind = toks.size() - 9 ;
                    nind *= numhap;
                    for (int i = 9 ; i < toks.size(); i++){
                        string id = toks[i];
                        indiv.push_back (ind(id,'U',"G"));
                        if (!isgenotype)
                            indiv.push_back (ind(id,'U',"G"));

                    }
                    setind = true;
                }
            }
        } else { 
            vector<string> toks;
            functions::tokenize (line.c_str(), toks, " \t");
            string id = toks[2];
            string chr = toks[0];
            double physpos =  atof(toks[1].c_str());
            double genpos = physpos*1.3/1e8;
            if (id.compare (".")==0) { 
                id = chr+":"+toks[1];
            }
            string ref = toks[3];
            string var = toks[4];
            string pass = toks[6];
            string desc = toks[7];
            string format = toks[8];
            if (pass.compare("PASS")!=0)
                continue;
            if (ref.compare (".")==0)
                continue;
            if (desc.compare (".") == 0 ) {
                cout << "Warning: No INFO field found. Assuming site is a SNP " << endl;
            }
            vector<string> tok1;
            functions::tokenize (desc.c_str(), tok1, ";");
            if (tok1[0].compare("SNP")!=0)
                continue;


            // Fill the snps
            if (chr.compare(prevchr)==0) {
                if (physpos >= prevppos){
                    prevgpos = genpos;
                    prevppos = physpos;
                } else {
                    cerr << "Invalid snps in  " << genofile <<" : Line " << linenum  << "\t Must have non-decreasing physical positions" << endl;
            //				exit (1);
                }
            } else {
                prevchr = chr;
                prevgpos = genpos;
                prevppos = physpos;
                if (snps.find(chr) != snps.end ()) {
                    cerr << "Invalid snps in "<< genofile << " : Line " << linenum << "\t SNPs on a chromosome must be consecutive" << endl;
                    exit (1);
                }
                chrs.push_back(chr);
                snps[chr] = vector<snp>();
//                snps[chr].reserve (size);
            }

/*            if (j==size){ 
                size += inc;
                snps[chr].reserve(size);
            }
            */

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

            // Fill in the genotypes
            tok1.resize (0);
            functions::tokenize (format.c_str(), tok1,":");
            int gind = -1;

            for (int i = 0 ; i < tok1.size(); i++){
                if (tok1[i].compare ("GT")==0) { 
                    gind = i;
                    break;
                }
            }
            if (gind == -1)
                continue;

            snps[chr].push_back (snp(id,chr,physpos,genpos,var,ref));
            vector<int> &gtype = snps[chr].back().gtype;
            j++;

            if (setind) {
                if (nind != numhap*(toks.size()-9)){
					cerr << "Bad genotype in " << genofile << ": Line "<< linenum  << endl;
					exit(1);
				}
            } else {
                nind = toks.size() - 9;
                nind *= numhap;
                setind = true;
            }

            gtype.resize(nind, 9);
            
            for (int i = 9, k = 0; i < toks.size(); i++, k++){ 
                tok1.resize(0);
                functions::tokenize(toks[i], tok1,",");
                if (isgenotype)  {
                    int a1 =  tok1[gind][0];
                    int a2 =  tok1[gind][2];
                    gtype[k] = (a1-48) + (a2-48);
                } else { 
                    int a1 =  tok1[gind][0];
                    int a2 =  tok1[gind][2];
                    gtype[2*k] = a1-48;
                    gtype[2*k+1] = a2-48;
                }
            }
            
        }
        if (j%10000==0) {
            io::println ("Read " + tostring (j) + " lines" ,0);
        }
    }
    nchr = chrs.size ();
    nsnps = j;
    
    if (io::debug >= 2 ) {
        cout << "genotype = " << to_string ();
    }

    ig.close ();

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
		io::println ("Length = "  + tostring(line.length()),4);
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
			if ( io::debug >= 4) {
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

	if (io::debug >= 3) { 
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
        indivmap[id] = (linenum-1);
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


/*
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

//	for (int i = 0 ; i < chrs.size(); i++){
//		cout << i << "\t" <<  chrs[i] << "\t" << snps[chrs[i]].size() << endl;
//	}

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
	
}*/

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

// Gets the index of the snp (given id)
int genotype::get_index (string id) { 
	if (allsnpsmap.find (id) == allsnpsmap.end())
		return -1;
	else {
		pair<string, int> p = allsnpsmap[id];
		return p.second;
	}
}

// Gets the chr, index of snp
pair<string,int> genotype::get_snp (string id) { 
	pair<string,int>  p ("",-1);
	if (allsnpsmap.find (id) != allsnpsmap.end())
		p = allsnpsmap[id];
	return p;
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
//				exit (1);
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

    int cumnum = 0;
    cumnumsnps.resize (nchr);
	for (int i  = 0 ; i < nchr; i++) { 
		numsnps[chrs[i]] = snps[chrs[i]].size();
        cumnumsnps[i] = cumnum;
        cumnum += numsnps[chrs[i]];
        
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
