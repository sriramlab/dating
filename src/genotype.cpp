#include "intervalmap.h"
#include "packedgtype.h"
#include "genotype.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "ind.h"
#include "gzstream.h"
#include "std.h"
#include "functions.h"
#include "fileio.h"


genotype::genotype (data *d) { 
	this->d = d;
    this->indivfilter = NULL;
    this->snpfilter =  NULL;
    this->prevsnpindex = -1;
	d->get_string ("genotypename",genofile,true);
	d->get_string ("snpname",snpfile,true);
	d->get_string ("indivname",indivfile,true);
	d->get_boolean ("isgenotype", isgenotype, true);
    goodsnpfile = "";
    givengoodsnp = d->get_string ("goodsnpname", goodsnpfile, false);
    givenbadsnp = d->get_string ("badsnpname", goodsnpfile, false);
    if (givengoodsnp && givenbadsnp) {
        cerr << "Can specify only one of goodsnpname and badsnpname"<<endl;
        exit(1);
    }

    indivfilterfile = "";
    givenindivfilter = d->get_string ("goodindname", indivfilterfile, false);
    d->get_boolean ("pack", pack, true);

    string inputformat = "eigenstrat";
    d->get_string ("inputformat", inputformat, false);
    d->get_boolean ("hasmissing", hasmissing, true);
    
    if (inputformat.compare ("prob")==0 && (!pack || hasmissing)){
        cerr << "Incompatible options for prob entry in genotype::genotype" << endl;
        cerr << "Need pack to be true"<<endl;
        cerr << "Need hasmissing to false" <<endl;
        exit(1);
    }
	
	read_snps (snpfile, goodsnpfile);
	read_ind (indivfile, indivfilterfile);
    if (pack) {
        pg = new packedgtype (nind, nsnp, cumnumsnps, isgenotype, hasmissing, NULL, NULL, inputformat);
        if ( inputformat.compare("eigenstrat")==0) 
            pg->read_es (genofile, indivfilter, snpfilter);
        else if (inputformat.compare("packedeigenstrat")==0)
            pg->read_packedes (genofile, indivfilter, snpfilter);
        else if (inputformat.compare("prob")==0)
            pg->read_es (genofile, indivfilter, snpfilter);
        else  {
            cerr << "Unknown format" << endl;
            exit(1);
        }
        nind = pg->nind ;
        set_freq ();
    } else  {
        read_es (genofile);
    }

	givenweights = d->get_string ("weightname",weightfile,false);
	if (givenweights) 
		read_weights (weightfile);
	d->get_boolean ("outputfail", outputfail, false, false);


}

genotype::genotype (string snpfile, string genofile, bool isgenotype, bool hasmissing, bool pack, string sep) { 
    this->pack = pack;
    this->sep = sep;
	read_snps (snpfile, "");
	this->isgenotype = isgenotype;
    this->indivfilter = NULL;
    this->snpfilter =  NULL;
    this->prevsnpindex = -1;
	givenind = false;

    if (pack) {
        pg = new packedgtype (nsnp, cumnumsnps, isgenotype, hasmissing);
        pg->read_es (genofile, NULL, NULL);
        set_freq ();
        nind = pg->nind ;
    } else  {
        read_es (genofile);
    }

}

genotype::genotype (string snpfile, string indivfile, string genofile, string format , bool isgenotype , string goodindname, string goodsnpname, bool hasmissing, bool pack, string sep) { 
    this->pack = pack;
    this->sep =sep;
    this->indivfilter = NULL;
    this->snpfilter =  NULL;
    this->isgenotype = isgenotype;
    this->prevsnpindex = -1;

    struct timeval now;
    gettimeofday(&now, NULL);
    long starttime = now.tv_sec * UMILLION + now.tv_usec;
    if (format.compare ("eigenstrat")==0) {
    	read_snps (snpfile, goodsnpname);

        if (functions::fileExists(indivfile)) {
            read_ind (indivfile, goodindname);
            
            if (pack) {
                pg = new packedgtype (nind, nsnp, cumnumsnps, isgenotype, hasmissing);
                pg->read_es (genofile, indivfilter, snpfilter);
                set_freq();
                
            } else  {
                read_es (genofile);
            }
        } else {
            givenind = false;
            read_es (genofile);
        }
    } else if (format.compare ("packedeigenstrat")==0){
        if (goodsnpname.compare("")==0 && goodindname.compare ("")==0) { 

            if (functions::fileExists (snpfile)){
                read_snps (snpfile, goodsnpname);
            } else {
                cerr << "Cannot find snpfile " << endl;
                exit (1);
            }
            if (functions::fileExists(indivfile)) {
                read_ind (indivfile, goodindname);
            } else {
                cerr << "Cannot find snpfile " << endl;
                exit (1);
            }


            pg = new packedgtype (nind, nsnp, cumnumsnps, isgenotype, hasmissing);
            pg->read_packedes (genofile);
            set_freq ();

        } else {

            if (functions::fileExists (snpfile)){
                read_snps (snpfile, goodsnpname);
            } else {
                cerr << "Cannot find snpfile " << endl;
                exit (1);
            }
            if (functions::fileExists(indivfile)) {
                read_ind (indivfile, goodindname);
            } else {
                cerr << "Cannot find snpfile " << endl;
                exit (1);
            }


            pg = new packedgtype (nind, nsnp, cumnumsnps, isgenotype, hasmissing);
            pg->read_packedes (genofile, indivfilter, snpfilter);
            set_freq ();
        }
    } else if (format.compare ("vcf")==0) { 
//        read_vcf (genofile, isgenotype, filterpolicy);
    }
    gettimeofday(&now, NULL);
    long endtime = now.tv_sec * UMILLION + now.tv_usec;
    long elapsed = endtime - starttime;
    elapsed/=1e6;
    if (io::debug >= 1)
        cout << "# Time = " << elapsed << endl;

}

void genotype::read_vcf (string genofile, bool isgenoytpe, intervalmap *imap, string goodsnpsname, bool isx) { 
    igzstream ig (genofile.c_str());
    string line;
	vector<snp> tmp;
	int j = 0 ;
	maxgpos = unordered_map<string, double> ();
	mingpos = unordered_map<string, double>();
	maxppos = unordered_map<string, double> ();
	minppos = unordered_map<string, double>();
    unordered_map <string, string> formatmap;
    unordered_map <string, string> infomap;
    unordered_map <string, double> attributemap;
	chrs = vector<string> ();
	double prevgpos ;
	double prevppos;
	string prevchr="";
    int size =  10000;
    int inc = 10000;
    
    unordered_map <string, string> * gsmap  = NULL;
    bool givengoodsnps  = false;
    vector<ind> tmpindiv;

    if (functions::fileExists (goodsnpsname) ) { 
        gsmap =new unordered_map <string, string> ();
        fileio::read_map (goodsnpsname, gsmap, 0);
        givengoodsnps  = true;
        cout << "goodsnpname =  " << goodsnpsname  << endl;
    }


    int linenum = 0 ;
    bool setind = false;
    int numhap = 2;
    int index;
    int intervalind;
    vector<interval> *intervalsonchr;
    if (isgenotype)
        numhap = 1;
    while (std::getline (ig, line) ) {
        linenum ++;
        char c = line[0];
        if (c=='#') { 
            if (line.find ("CHROM")!=string::npos) {
                if (!isx && !givenind)  {
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
                } else if (isx) { 
                    vector<string> toks;
                    functions::tokenize (line.c_str(), toks, " \t");
                    for (int i = 9 ; i < toks.size(); i++){
                        string id = toks[i];
                        tmpindiv.push_back (ind(id,'U',"G"));
                    }
                }
            } else { 
                if (line.find ("INFO=")!=string::npos) { 
                    vector<string> toks;
                    functions::tokenize (line.c_str(), toks, ",");
                    string key;
                    string val;
                    for (int i = 0 ; i < toks.size (); i++)  {
                        vector<string> toks1;
                        functions::tokenize (toks[i].c_str(), toks1, "=");
                        if (toks[i].find ("ID=")!=string::npos) {
                            key = toks1[1];
                        } else if (toks[i].find ("Type=") !=string::npos) { 
                            val = toks1[1];
                        }
                    }
                    infomap[key] = val;
                } else if (line.find ("FORMAT=")!=string::npos) { 
                    vector<string> toks;
                    functions::tokenize (line.c_str(), toks, ",");
                    string key;
                    string val;
                    for (int i = 0 ; i < toks.size (); i++)  {
                        vector<string> toks1;
                        functions::tokenize (toks[i].c_str(), toks1, "=");
                        if (toks[i].find ("ID=")!=string::npos) {
                            key = toks1[1];
                        } else if (toks[i].find ("Type=") !=string::npos) { 
                            val = toks1[1];
                        }
                    }
                    formatmap[key] = val;
                }
            }
        } else { 
            vector<string> toks;
            functions::tokenize (line.c_str(), toks, " \t");
            string id = toks[2];
            string chr = toks[0];
            double physpos =  atof(toks[1].c_str());
            double genpos = physpos*1.3/1e8;

            if ( chr.compare(prevchr)!=0 ) {
                intervalsonchr = &(imap->intervals[chr]);
                intervalind = 0 ;
                for (int i = 0 ; i < intervalsonchr->size() ; i++) { 
                    interval &tmpint = (*intervalsonchr)[i];
                    if (tmpint.end > physpos) {
                        intervalind = i;
                        break;
                    }
                }
            } else {
                for (int i = intervalind; i < intervalsonchr->size(); i++)  {
                    interval &tmpint = (*intervalsonchr)[i];
                    if (tmpint.end > physpos) {
                        intervalind = i;
                        break;
                    }
                }
            }
            if (io::debug >= 0 ) {
                interval &tmpint = (*intervalsonchr)[intervalind];
                cout.precision (0); cout.setf(ios::fixed); 
                cout << noshowpoint <<  intervalind  << "\t" << tmpint.chr << "\t" << tmpint.start << "\t" << tmpint.end << endl;
            }

            if (isx) {
                // PAR
                if (physpos < 2699555 || physpos > 154930277)
                    continue;
            }


            if (id.compare (".")==0) { 
                id = chr+":"+toks[1];
            }

            id = chr + ":" + toks[1];
            string id1 = id;
            string ref = toks[3];
            string var = toks[4];
            string pass = toks[6];
            string info = toks[7];
            string format = toks[8];
            string geno = toks[9];

            vector<string> tok1;
            vector<string> tok2;
            vector<string> tok3;
            functions::tokenize (info.c_str(), tok1, ";");


            attributemap.clear ();
            for (int i =  0 ; i < tok1.size(); i++){ 
                tok2.resize(0);
                functions::tokenize (tok1[i], tok2, "=");
                if ( tok2.size() > 1)
                    attributemap[tok2[0]] = atof(tok2[1].c_str());
            }
            tok3.resize(0);tok2.resize(0);
            functions::tokenize (format.c_str(), tok3, ":");
            functions::tokenize (geno.c_str(), tok2, ":");
            for (int i  = 0 ; i < tok3.size(); i++) { 
                if (functions::isnumeric ( tok2[i].c_str())){
                    attributemap[tok3[i]] = atof(tok2[i].c_str());
                }
            }

            interval &tmpint = (*intervalsonchr)[intervalind];
            if (tmpint.start >physpos) {
                if (io::debug >= 1)
                    cout << "Dropping position " << id1 << ": not found in interval " <<endl;
                continue;
            }
            if (attributemap.find("GQ")!=attributemap.end()){
                if (attributemap["GQ"]>=30);
                else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id1 << ": below GQ cutoff"<<endl;
                    continue;
                }
            } else {
                if (io::debug >= 1)
                    cout << "Dropping position " << id1 << ": no GQ field"<<endl;
                continue;
            }

            if (ref.compare(".")==0) 
                ref="X";
            if (var.compare(".")==0)
                var="X";

            if (givengoodsnps) {
                if (gsmap->find(id1) == gsmap->end() || (*gsmap)[id1].compare ("1")==0 ) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id1 << ": not found in " << goodsnpsname<<endl;
                    continue;
                } else 
                    (*gsmap)[id1] = "1";
            }

            if (ref.compare (".")==0)
                continue;
            if (info.compare (".") == 0 ) {
                cout << "Warning: No INFO field found. Assuming site is a SNP " << endl;
            }

            cout << "Accepted " << id1 << "\t" << physpos << endl;



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

            snps[chr].push_back (snp(id,chr,physpos,genpos,var[0],ref[0]));
            vector<int> &gtype = snps[chr].back().gtype;
            j++;

            if (setind) {
                if (!isx && nind != numhap*(toks.size()-9)){
					cerr << "Bad genotype in " << genofile << ": Line "<< linenum  << endl;
					exit(1);
				}
            } else {
                if (!isx) {
                    nind = toks.size() - 9;
                    nind *= numhap;
                    setind = true;
                } else {
                    nind =  0 ;
                    for (int i = 9, k = 0; i < toks.size(); i++, k++){ 
                        tok1.resize(0);
                        functions::tokenize(toks[i], tok1,",");
                        nind ++;
                        indiv.push_back (tmpindiv[k]);
                        if (tok1[gind][1]=='|' && !isgenotype) {
                            nind ++;
                            indiv.push_back (tmpindiv[k]);
                        }  
                    }
                    setind = true;
                }

            }

            gtype.resize(nind, 9);
            
            for (int i = 9, k = 0, l = 0 ; i < toks.size(); i++, k++){ 
                tok1.resize(0);
                functions::tokenize(toks[i], tok1,",");
                if (!isx) { 
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
                } else { 
                    if (tok1[gind][1]=='|') {
                        if (isgenotype ) { 
                            int a1 =  tok1[gind][0];
                            int a2 =  tok1[gind][2];
                            gtype[k] = (a1-48) + (a2-48);
                        } else { 
                            int a1 =  tok1[gind][0];
                            int a2 =  tok1[gind][2];
                            gtype[l++] = a1-48;
                            gtype[l++] = a2-48;
                        }
                    } else {
                        if (isgenotype ) { 
                            int a1 =  tok1[gind][0];
                            gtype[k] = a1-48;
                        } else { 
                            int a1 =  tok1[gind][0];
                            gtype[l++] = a1-48;
                        }
                    }
                }
            }
            
        }
        if (j%10000==0) {
            io::println ("Read " + tostring (j) + " lines" ,0);
        }
    }
    nchr = chrs.size ();
    nsnp = j;
    
    if (io::debug >= 2 ) {
        cout << "genotype = " << to_string ();
    }

    ig.close ();

}

void genotype::read_vcf (string genofile, bool isgenoytpe, string filterpolicy, string goodsnpsname, bool isx) { 
    igzstream ig (genofile.c_str());
    string line;
	vector<snp> tmp;
	int j = 0 ;
	maxgpos = unordered_map<string, double> ();
	mingpos = unordered_map<string, double>();
	maxppos = unordered_map<string, double> ();
	minppos = unordered_map<string, double>();
    unordered_map <string, string> formatmap;
    unordered_map <string, string> infomap;
    unordered_map <string, double> attributemap;
	chrs = vector<string> ();
	double prevgpos ;
	double prevppos;
	string prevchr="";
    int size =  10000;
    int inc = 10000;
    
    unordered_map <string, string> * gsmap  = NULL;
    bool givengoodsnps  = false;
    vector<ind> tmpindiv;

    if (functions::fileExists (goodsnpsname) ) { 
        gsmap =new unordered_map <string, string> ();
        fileio::read_map (goodsnpsname, gsmap, 0);
        givengoodsnps  = true;
        cout << "goodsnpname =  " << goodsnpsname  << endl;
    }


    int linenum = 0 ;
    bool setind = false;
    int numhap = 2;
    int index;
    int intervalind;
    vector<interval> *intervalsonchr;
    if (isgenotype)
        numhap = 1;
    while (std::getline (ig, line) ) {
        linenum ++;
        char c = line[0];
        if (c=='#') { 
            if (line.find ("CHROM")!=string::npos) {
                if (!isx && !givenind)  {
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
                } else if (isx) { 
                    vector<string> toks;
                    functions::tokenize (line.c_str(), toks, " \t");
                    for (int i = 9 ; i < toks.size(); i++){
                        string id = toks[i];
                        tmpindiv.push_back (ind(id,'U',"G"));
                    }
                }
            } else { 
                if (line.find ("INFO=")!=string::npos) { 
                    vector<string> toks;
                    functions::tokenize (line.c_str(), toks, ",");
                    string key;
                    string val;
                    for (int i = 0 ; i < toks.size (); i++)  {
                        vector<string> toks1;
                        functions::tokenize (toks[i].c_str(), toks1, "=");
                        if (toks[i].find ("ID=")!=string::npos) {
                            key = toks1[1];
                        } else if (toks[i].find ("Type=") !=string::npos) { 
                            val = toks1[1];
                        }
                    }
                    infomap[key] = val;
                } else if (line.find ("FORMAT=")!=string::npos) { 
                    vector<string> toks;
                    functions::tokenize (line.c_str(), toks, ",");
                    string key;
                    string val;
                    for (int i = 0 ; i < toks.size (); i++)  {
                        vector<string> toks1;
                        functions::tokenize (toks[i].c_str(), toks1, "=");
                        if (toks[i].find ("ID=")!=string::npos) {
                            key = toks1[1];
                        } else if (toks[i].find ("Type=") !=string::npos) { 
                            val = toks1[1];
                        }
                    }
                    formatmap[key] = val;
                }
            }
        } else { 
            vector<string> toks;
            functions::tokenize (line.c_str(), toks, " \t");
            string id = toks[2];
            string chr = toks[0];
            double physpos =  atof(toks[1].c_str());
            double genpos = physpos*1.3/1e8;


            if (isx) {
                // PAR
                if (physpos < 2699555 || physpos > 154930277)
                    continue;
            }


            if (id.compare (".")==0) { 
                id = chr+":"+toks[1];
            }
            id = chr + ":" + toks[1];

            string ref = toks[3];
            string var = toks[4];
            string pass = toks[6];
            string info = toks[7];
            string format = toks[8];
            string geno = toks[9];



            vector<string> tok1;
            vector<string> tok2;
            vector<string> tok3;
            functions::tokenize (info.c_str(), tok1, ";");


            attributemap.clear ();
            for (int i =  0 ; i < tok1.size(); i++){ 
                tok2.resize(0);
                functions::tokenize (tok1[i], tok2, "=");
                if ( tok2.size() > 1)
                    attributemap[tok2[0]] = atof(tok2[1].c_str());
            }
            tok3.resize(0);tok2.resize(0);
            functions::tokenize (format.c_str(), tok3, ":");
            functions::tokenize (geno.c_str(), tok2, ":");
            for (int i  = 0 ; i < tok3.size(); i++) { 
                if (functions::isnumeric ( tok2[i].c_str())){
                    attributemap[tok3[i]] = atof(tok2[i].c_str());
                }
            }



            if (filterpolicy.compare ("160812")==0) { 
                id = chr+":"+toks[1];

                if (ref.length() > 1 || var.length() > 1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": indel" <<endl;
                    continue;
                } 

                if (attributemap.find("Map20")!=attributemap.end()){
                    if (attributemap["Map20"] !=1) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": below Map20 cutoff"<<endl;
                        continue;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no Map20 field"<<endl;
                    continue;
                }
                if (attributemap.find("GQ")!=attributemap.end()){
                    if (attributemap["GQ"]>=30);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": below GQ cutoff"<<endl;
                        continue;
                    }
                } else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no GQ field"<<endl;
                    continue;
                }
                if (attributemap.find ("DP")!=attributemap.end()){
                    double rd = attributemap["DP"];
                    if (rd >= 31 && rd <= 79);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": outside DP range"<<endl;
                        continue;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no DP field"<<endl;
                    continue;
                }

                if (ref.compare(".")==0) 
                    ref="X";
                if (var.compare(".")==0)
                    var="X";
            
                if (givengoodsnps) {
                    if (gsmap->find(id) == gsmap->end() || (*gsmap)[id].compare ("1")==0 ) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": not found in " << goodsnpsname<<endl;
                        continue;
                    } else 
                        (*gsmap)[id] = "1";
                }
            } else if (filterpolicy.compare ("010413")==0) { 
                id = chr+":"+toks[1];

                if (ref.length() > 1 || var.length() > 1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": indel" <<endl;
                    continue;
                } 

                if (attributemap.find("GQ")!=attributemap.end()){
                    if (attributemap["GQ"]>=30);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": below GQ cutoff"<<endl;
                        continue;
                    }
                } else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no GQ field"<<endl;
                    continue;
                }
                if (attributemap.find ("DP")!=attributemap.end()){
                    double rd = attributemap["DP"];
                    if (rd >= 31 && rd <= 79);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": outside DP range"<<endl;
                        continue;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no DP field"<<endl;
                    continue;
                }

                if (ref.compare(".")==0) 
                    ref="X";
                if (var.compare(".")==0)
                    var="X";
            
                if (givengoodsnps) {
                    if (gsmap->find(id) == gsmap->end() || (*gsmap)[id].compare ("1")==0 ) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": not found in " << goodsnpsname<<endl;
                        continue;
                    } else 
                        (*gsmap)[id] = "1";
                }
            } else if (filterpolicy.compare ("1kg")==0 ) {
                if (pass.compare ("PASS")!=0)
                    continue;
                if (tok1[0].compare("SNP")!=0)
                    continue;
            } else if (filterpolicy.compare ("gonl")==0) {
                if (pass.compare ("PASS")!=0)
                    continue;
                id = chr+":"+toks[1];
                if (ref.length() > 1 || var.length() > 1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": indel" <<endl;
                    continue;
                }

                if (ref.compare(".")==0) 
                    ref="X";
                if (var.compare(".")==0)
                    var="X";
            
                if (givengoodsnps) {
                    if (gsmap->find(id) == gsmap->end() || (*gsmap)[id].compare ("1")==0 ) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": not found in " << goodsnpsname<<endl;
                        continue;
                    } else 
                        (*gsmap)[id] = "1";
                }

            } else if (filterpolicy.compare ("051012")==0) { 
                id = chr+":"+toks[1];

                if (ref.length() > 1 || var.length() > 1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": indel" <<endl;
                    continue;
                } 

                if (attributemap.find("Map20")!=attributemap.end()){
                    if (attributemap["Map20"] !=1) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": below Map20 cutoff"<<endl;
                        continue;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no Map20 field"<<endl;
                    continue;
                }
                if (attributemap.find("GQ")!=attributemap.end()){
                    if (attributemap["GQ"]>=30);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": below GQ cutoff"<<endl;
                        continue;
                    }
                } else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no GQ field"<<endl;
                    continue;
                }
                if (attributemap.find ("DP")!=attributemap.end()){
                    double rd = attributemap["DP"];
                    if (rd >= 16 && rd <= 46);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": outside DP range"<<endl;
                        continue;
                    }
                } else  {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no DP field"<<endl;
                    continue;
                }

                if (ref.compare(".")==0) 
                    ref="X";
                if (var.compare(".")==0)
                    var="X";
            
                if (givengoodsnps) {
                    if (gsmap->find(id) == gsmap->end() || (*gsmap)[id].compare ("1")==0 ) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": not found in " << goodsnpsname<<endl;
                        continue;
                    } else 
                        (*gsmap)[id] = "1";
                }
            } else { 

                if (ref.length() > 1 || var.length() > 1) {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": indel" <<endl;
                    continue;
                }

                if (attributemap.find("GQ")!=attributemap.end()){
                    if (attributemap["GQ"]>=30);
                    else {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": below GQ cutoff"<<endl;
                        continue;
                    }
                } else {
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": no GQ field"<<endl;
                    continue;
                }

                if (ref.compare(".")==0) 
                    ref="X";
                if (var.compare(".")==0)
                    var="X";
 
                if (givengoodsnps) {
                    if (gsmap->find(id) == gsmap->end() || (*gsmap)[id].compare ("1")==0 ) {
                        if (io::debug >= 1)
                            cout << "Dropping position " << id << ": not found in " << goodsnpsname<<endl;
                        continue;
                    } else 
                        (*gsmap)[id] = "1";
                }
            }

            if (ref.compare (".")==0)
                continue;
            if (info.compare (".") == 0 ) {
                cout << "Warning: No INFO field found. Assuming site is a SNP " << endl;
            }



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
            }


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

            snps[chr].push_back (snp(id,chr,physpos,genpos,var[0],ref[0]));
            vector<int> &gtype = snps[chr].back().gtype;
            j++;

            if (setind) {
                if (!isx && nind != numhap*(toks.size()-9)){
					cerr << "Bad genotype in " << genofile << ": Line "<< linenum  << endl;
					exit(1);
				}
            } else {
                if (!isx) {
                    nind = toks.size() - 9;
                    nind *= numhap;
                    setind = true;
                } else {
                    nind =  0 ;
                    for (int i = 9, k = 0; i < toks.size(); i++, k++){ 
                        tok1.resize(0);
                        functions::tokenize(toks[i], tok1,",");
                        nind ++;
                        indiv.push_back (tmpindiv[k]);
                        if (tok1[gind][1]=='|' && !isgenotype) {
                            nind ++;
                            indiv.push_back (tmpindiv[k]);
                        }  
                    }
                    setind = true;
                }

            }

            gtype.resize(nind, 9);
            
            for (int i = 9, k = 0, l = 0 ; i < toks.size(); i++, k++){ 
                tok1.resize(0);
                functions::tokenize(toks[i], tok1,",");
                if (!isx) { 
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
                } else { 
                    if (tok1[gind][1]=='|') {
                        if (isgenotype ) { 
                            int a1 =  tok1[gind][0];
                            int a2 =  tok1[gind][2];
                            gtype[k] = (a1-48) + (a2-48);
                        } else { 
                            int a1 =  tok1[gind][0];
                            int a2 =  tok1[gind][2];
                            gtype[l++] = a1-48;
                            gtype[l++] = a2-48;
                        }
                    } else {
                        if (isgenotype ) { 
                            int a1 =  tok1[gind][0];
                            gtype[k] = a1-48;
                        } else { 
                            int a1 =  tok1[gind][0];
                            gtype[l++] = a1-48;
                        }
                    }
                }
            }
            
        }
        if (j%10000==0) {
            io::println ("Read " + tostring (j) + " lines" ,0);
        }
    }
    nchr = chrs.size ();
    nsnp = j;
    
    if (io::debug >= 2 ) {
        cout << "genotype = " << to_string ();
    }

    ig.close ();

}


void genotype::set_weights (string weightfile){ 
	givenweights  = true;
	read_weights (weightfile);
}

void genotype::write_genotypes(string filename) {
    if (pack)
        pg->write_packedes (filename);
}

ostream& genotype::write_genotypes ( ostream & os) {
    if (pack) 
        return pg->write_packedes (os);
}

/*
ostream& genotype::write_genotypes ( ostream & os, vector<bool> *snpfilter) {
    if (pack) 
        return pg->write_packedes (os, snpfilter);
}*/



void genotype::read_packedes (string filename)  {
    pg->read_packedes (filename);

    // Run checks
    if (pg->cumnumsnps.size() != nchr){
        cerr << "Number of chromosomes do not match in" << filename << endl;
        exit (1);
    }
    if (pg->nind != nind) { 
        cerr << "Number of individuals do not match in" << filename << endl;
        exit (1);
    }
    if (pg->nsnp != nsnp) { 
        cerr << "Number of SNPs do not match in" << filename << endl;
        exit (1);
    }
}

void genotype::read_es (string filename) {
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	int chrindex = 0;
	int snpindex = 0;
    int effectiveindex  = 0 ;
	bool setind = false;
	while ( std::getline (inp, line)){
		linenum ++;
        if (linenum > originalnsnp ) {
            cerr << "Genotype file " << filename <<  " has more SNPs than the SNP file" << endl;
            exit (1);
        }
		char c = line[0];
		if (c=='#')
			continue;
		
		if (line.empty())
			continue;

        effectiveindex ++;
        if (snpfilter != NULL && (*snpfilter)[effectiveindex-1]==false) 
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
			if (originalnind != line.length()) {
				cerr << "Bad genotype in " << filename << ": Line "<< linenum  << endl;
				exit(1);
			}
		}

        unsigned int i =  0 ;
		for (int j = 0 ; j < line.length();j++){
            if (indivfilter!=NULL && (*indivfilter)[j]==false)
                continue;
			char c = line.at(j);
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
            i ++;
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

void genotype::read_ind (string filename, string filterfile) { 
	ifstream inp (filename.c_str());
    unordered_map <string, string> *tmpmap =  NULL;
    if (functions::fileExists (filterfile) ) { 
        tmpmap = new unordered_map<string,string>();
        fileio::read_map (filterfile,  tmpmap, 0);
    }

    indivfilter = new vector<bool> ();
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	vector<snp> tmp;
	int linenum  = 0;
    originalnind = 0;
	while ( std::getline (inp, line)){
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
        originalnind ++;
        if (tmpmap != NULL) {
            if (tmpmap->find(id)!=tmpmap->end()) {
                linenum ++;
                indiv.push_back (ind (id,g,group));
                indivmap[id] = (linenum-1);
                (*indivfilter).push_back(true);
                popindexmap[group].push_back (linenum-1);
            } else {
                (*indivfilter).push_back(false);
            }
        } else {
            linenum ++;
            indiv.push_back (ind (id,g,group));
            indivmap[id] = (linenum-1);
            (*indivfilter).push_back(true);
            popindexmap[group].push_back (linenum-1);
        }
    }
	nind  = indiv.size();
    prevgtype.resize(nind,0);

    if (nind == 0  ) {
        cerr << "No individuals after filtering" << endl;
        exit(1);
    }
	inp.close ();

	for (int i = 0 ; i < nchr; i++){
        vector<snp> &tmpsnps =  snps[chrs[i]];
        if (!pack) { 
            for (int  j = 0 ; j < tmpsnps.size(); j++){
                tmpsnps[j].gtype.resize(nind,9);
            }
        }
	}
	givenind = true;
    delete tmpmap;
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


void genotype::read_snps (string filename, string goodsnpsname){
    snpfilter = new vector <bool> ();
    unordered_map <string, string> * gsmap  = NULL;
    bool givensnpstofilter  = false;
    vector<ind> tmpindiv;

    if (functions::fileExists (goodsnpsname) ) { 
        gsmap =new unordered_map <string, string> ();
        fileio::read_map (goodsnpsname, gsmap, 0);
        givensnpstofilter  = true;
    } 

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
    originalnsnp = 0 ;
	
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

        (*snpfilter).push_back (true);
        originalnsnp ++;

        if (givensnpstofilter) {
            if (givengoodsnp) {
                if (gsmap->find(id) == gsmap->end() ) {
                    (*snpfilter)[originalnsnp-1] = false;
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": not found in " << goodsnpsname<<endl;
                    continue;
                }
            } else {
                if (gsmap->find(id) != gsmap->end()  ) {
                    (*snpfilter)[originalnsnp-1] = false;
                    if (io::debug >= 1)
                        cout << "Dropping position " << id << ": found in " << goodsnpsname<<endl;
                    continue;
                }

            }
        }


		if (chr.compare(prevchr)==0) {
			if (physpos >= prevppos){
				prevgpos = genpos;
				prevppos = physpos;
			} else {
				cerr << "Invalid snp file " << filename <<" : Line " << linenum  << "\t Must have non-decreasing physical positions" << endl;
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
	nsnp = j;
	nchr = chrs.size();

    if (nsnp == 0 ) {
        cerr << "No snps after filtering" << endl;
        exit(1);
    }

    int cumnum = 0;
    cumnumsnps.resize (nchr);
	for (int i  = 0 ; i < nchr; i++) { 
        chrindexmap[chrs[i]] = i;
		numsnps[chrs[i]] = snps[chrs[i]].size();
        cumnumsnps[i] = cumnum;
        cumnum += numsnps[chrs[i]];
        
        if (!pack)
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
	nsnp = j;
	nchr = chrs.size();

    int cumnum = 0;
    cumnumsnps.resize (nchr);
	for (int i  = 0 ; i < nchr; i++) { 
		numsnps[chrs[i]] = snps[chrs[i]].size();
        cumnumsnps[i] = cumnum;
        cumnum += numsnps[chrs[i]];
        
        if (!pack)
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


void genotype::set_freq ( )  {
    if (!pack || pg==NULL) { 
        cerr << "Erroneous call to set_freq" << endl;
        exit (1);
    }
    for (int i =  0 ; i < chrs.size(); i++) {
        vector<snp> &s = snps[chrs[i]];
        pg->set_freq (i, s);
    }
}


double genotype::get_freq (int chrindex, unsigned int snpindex, string popname) {
    if (snpindex != prevsnpindex) { 
        get_geno (chrindex, snpindex, prevgtype);
        prevsnpindex = snpindex;
    }
    return get_freq (prevgtype, popname);

}

double genotype::get_freq (unsigned int snpindex, string popname) {
    if (snpindex != prevsnpindex) { 
        get_geno (snpindex, prevgtype);
        prevsnpindex = snpindex;
    }
    return get_freq (prevgtype, popname);
}

double genotype::get_freq (unsigned int snpindex, int &vc, int &rc, int &mc, string popname) { 
    if (snpindex != prevsnpindex) { 
        get_geno (snpindex, prevgtype);
        prevsnpindex = snpindex;
    }
    return get_freq (prevgtype, vc, rc, mc, popname);
}

double genotype::get_freq (int chrindex, unsigned int snpindex, int &vc, int &rc, int &mc, string popname) { 
    if (snpindex != prevsnpindex) { 
        get_geno (chrindex, snpindex, prevgtype);
        prevsnpindex = snpindex;
    }
    return get_freq (prevgtype, vc, rc, mc, popname);
}



double genotype::get_freq (vector<int> &gtype, string popname) { 
    double f  = -1;
    int fnum = 0; int fdenom =  0;
    if (popindexmap.find(popname)!=popindexmap.end()){
        vector<int> &ind = popindexmap[popname];
        for (int i  = 0 ;  i < ind.size();i++){
            int g = gtype[ind[i]];
            if (g==9)
                continue;
            fnum += g;
            fdenom += isgenotype?2:1;
        } 
        f = fdenom>0?(1.*fnum)/fdenom:-1;
        return f;
    } else if (popname.compare("")==0){
        for (int i  = 0 ;  i < nind;i++){
            int g = gtype[i];
            if (g==9)
                continue;
            fnum += g;
            fdenom += isgenotype?2:1;
        } 
        f = fdenom>0?(1.*fnum)/fdenom:-1;
        return f;
    } else {
        cerr << "Cannot find population " << popname << endl;
        exit(1);
    }
}


double genotype::get_freq (vector<int> &gtype, int &vc, int &rc, int &mc, string popname) { 
    double f  = -1;
    int fnum = 0; int fdenom =  0;
    if (popindexmap.find(popname)!=popindexmap.end()){ 
        vector<int> &ind = popindexmap[popname];
        vc = 0; rc = 0 ; mc = 0;
        int numalleles = isgenotype?2:1;

        for (int i  = 0 ;  i < ind.size();i++){
            int v = gtype[ind[i]];
            if (v==9) {
                mc += numalleles;
                continue;
            }
            int r = numalleles - v;
            vc += v;  rc += r;
        } 
        fdenom = vc + rc; fnum = vc;
        f = fdenom>0?(1.*fnum)/fdenom:-1;
        return f;
    } else if (popname.compare ("")==0){
        vc = 0; rc = 0 ; mc = 0;
        int numalleles = isgenotype?2:1;

        for (int i  = 0 ;  i < nind; i++){
            int v = gtype[i];
            if (v==9) {
                mc += numalleles;
                continue;
            }
            int r = numalleles - v;
            vc += v;  rc += r;
        } 
        fdenom = vc + rc; fnum = vc;
        f = fdenom>0?(1.*fnum)/fdenom:-1;
        return f;
    } else {
        cerr << "Cannot find population " << popname << endl;
        exit(1);
    }
}


