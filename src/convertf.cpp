#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "convertf.h"

ofstream outfs;

convertf::convertf (int argc, char *argv[]) {
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
    string inputformat = "eigenstrat";
	p->get_string ("genotypename",geno1,true);
	p->get_string ("snpname",snp1,true);
	p->get_string ("indivname",ind1,true);
    vector<string> tok;
    functions::tokenize(geno1.c_str(), tok,".");
    if (tok[tok.size()-1].compare("gz")==0){
        if (tok[tok.size()-2].compare("vcf")==0) {
            inputformat = "vcf";
            filterpolicy = "1kg";
            p->get_string ("filterpolicy", filterpolicy, false);
        }
    }
    cout << "Inputformat = " << inputformat << endl;

    outputformat = "eigenstrat";
    p->get_string ("outputformat", outputformat, false);
	p->get_boolean ("isgenotype", isgenotype, true);
	p->get_boolean ("isx", isx, false);
    p->get_int ("debug" , io::debug, 0 );


    extracttrios = false;
    trioflag = true;
    p->get_boolean ("extracttrios", extracttrios, false);

    string gsfile =  "";
	bool goodsnps = p->get_string ("goodsnpname",gsfile,false);

    if (io::debug>=2) 
        cout << "Reading genotypes from " <<geno1<<endl;
    g = new genotype (snp1,ind1, geno1, inputformat, isgenotype, "", "", false, false, "");

    if (io::debug>=2) 
        cout << "Read genotypes from " <<geno1<<endl;

    if (inputformat.compare ("vcf") ==0 ) {
        intervalmap *imap;
        unordered_map <string, vector<interval> > *intervals;
        string chr;
        bool givenchr = p->get_string ("chr", chr, false );
        if (filterpolicy.compare ("160812")!=0 && filterpolicy.compare ("010413")!=0
                &&filterpolicy.compare ("1kg")!=0 &&filterpolicy.compare ("051012")!=0 
                && filterpolicy.compare ("basic")!=0 && filterpolicy.compare ("gonl")!=0 ) { 
            imap = new intervalmap (filterpolicy, chr);

            if (io::debug >= 2){ 
                cout << "Read interval file " << filterpolicy << endl;

                vector<interval>  &intervalsonchr = imap->intervals["22"];
                for (int i = 0 ; i < intervalsonchr.size (); i++) {
                    interval &tmpint = intervalsonchr[i];
                    cout << tmpint.chr << "\t" << tmpint.start <<"\t" << tmpint.end << endl;

                }
            }

            g->read_vcf (geno1, isgenotype, imap , gsfile, isx);
        } else 
            g->read_vcf (geno1, isgenotype, filterpolicy , gsfile, isx);
    }

}

void convertf::convert ()  {

	p->get_string("snpoutname|snpoutfilename",snpoutfile,true);
	p->get_string("indivoutname|indoutfilename",indivoutfile,true);
	p->get_string("genotypeoutname|genooutfilename",genotypeoutfile,true);

	string gsfile;
	string flipfile;
	string polarizefile;
    bool useonlyconfidentsites;
	string filterfile;
	string newindivnamefile;
	string chr;
    string gstart;
    string gend;
    string gstartchr;
    double gstartpos;
    string gendchr;
    double gendpos;


	bool goodsnps = p->get_string ("goodsnpname",gsfile,false);
    bool goodsnpstart = p->get_string ("goodsnpstart",gstart, false);
    if (goodsnpstart) {
        vector<string> tokens;
        functions::tokenize(gstart,tokens,": \t");
        gstartchr =  tokens[0];
        gstartpos = atof(tokens[1].c_str());
    }
    bool goodsnpend = p->get_string ("goodsnpend", gend, false);
    if (goodsnpend) {
        vector<string> tokens;
        functions::tokenize(gend,tokens,": \t");
        gendchr =  tokens[0];
        gendpos = atof(tokens[1].c_str());
    }

	bool badsnps = p->get_string ("badsnpname",gsfile,false);
	bool flipsnps = p->get_string ("flipsnpname", flipfile,false);
	bool polarize = p->get_string ("polarizesnpname", polarizefile, false);
    p->get_boolean("useonlyconfidentsites", useonlyconfidentsites, true, false); 
	bool filter = p->get_string ("filtersnpname", filterfile, false);
	bool pickchr = p->get_string ("chrom",chr,false);
	bool newindivname = p->get_string ("newindivname", newindivnamefile, false);
    bool keeporder  = true;
    p->get_boolean ("keeporder", keeporder, true, false);

	p->get_boolean ("uniqpos", uniqpos, false, false);

	if (goodsnps && badsnps) {
		cerr << "Can specify only one of goodsnpname or badsnpname" <<endl;
        exit(1);
	}

	snpmap *smap;
	unordered_map<string, string > *smapsnps = NULL;
	unordered_map<string, string > *pmapsnps = NULL;
	unordered_map<string, string > *filtermapsnps = NULL;
	unordered_map<string, string > *fmapsnps = NULL;
	unordered_map<string, string > *newindmap = NULL;
    vector<string> *newindvector = NULL;
//	unordered_set<snp,boost::hash<snp> > *smapsnps ;
	if (goodsnps ||badsnps){
		smapsnps = new unordered_map <string, string > ();
		fileio::read_map (gsfile, smapsnps,0);
		if (io::debug >= 2)  {
			typedef unordered_map<string, string>::iterator iter;
			for ( iter i  = smapsnps->begin();  i != smapsnps->end(); i++){
				cout << i->first << "\t" << i->second << endl;
			}
		}
	}

	if (polarize){  
		pmapsnps = new unordered_map <string, string> ();
		fileio::read_map (polarizefile,  pmapsnps, 0 , 1);
		if (io::debug >= 2)  {
			typedef unordered_map<string, string>::iterator iter;
			for ( iter i  = pmapsnps->begin();  i != pmapsnps->end(); i++){
				cout << i->first << "\t" << i->second << endl;
			}
		}
	}

	if (filter){  
		filtermapsnps = new unordered_map <string, string> ();
		fileio::read_map (filterfile, filtermapsnps, 0);
	}


	if (flipsnps){  
		fmapsnps = new unordered_map <string, string> ();
		fileio::read_map (flipfile, fmapsnps, 0 );
	}

	if (newindivname) { 
		newindmap = new unordered_map <string, string> ();
        newindvector = new vector<string> ();
		fileio::read_map (newindivnamefile, newindmap, 0 );
        if (!keeporder){
            fileio::read_vector(newindivnamefile, newindvector, 0);
        }
	}

	ofstream snpfs (snpoutfile.c_str());
	ofstream indivfs (indivoutfile.c_str());
	ofstream genofs (genotypeoutfile.c_str());


	int nind = g->nind;
	int nchr = g->nchr;
    vector< vector <int> > outputgeno;
    int outputsnps = 0 ;
    if (outputformat.compare ("ihs")==0) {
        if (isgenotype) {
            cerr << "Input genotypes need to be phased" << endl;
            exit(1);
        }
        for (int i = 0 ; i < nind; i++){
            outputgeno.push_back(vector<int> (g->nsnp));
        }
    }

	for (int i = 0 ; i < nchr; i++){
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];
        cout << "chr :"  << i << "\t" << tmpsnps.size()<< endl;
		if (pickchr && chr!=g->chrs[i])
			continue;
		double prevgpos ;
		bool flag = false;
		for (int j  = 0 ; j < tmpsnps.size();j++) {
			snp &s = tmpsnps[j];
            if (io::debug >= 2 ) {
                cout  << j << "\t" << s.id << endl;
            }
			if (goodsnps){
				if (smapsnps->find(s.id) == smapsnps->end())
					continue;
			} else if (badsnps)  {
				if (smapsnps->find(s.id) != smapsnps->end())
					continue;
			}
            if (goodsnpstart || goodsnpend) { 
                
                if (goodsnpstart) {
                    if (s.chr!=gstartchr || s.physpos < gstartpos)
                        continue;
                }
                if (goodsnpend) { 
                    if (s.chr!=gendchr || s.physpos > gendpos)
                        continue;
                }

            }

            if (polarize) {
                vector<int> &gtype = s.gtype;
                if ( pmapsnps->find (s.id) != pmapsnps->end ()) {
                    string a = (*pmapsnps)[s.id];
                    if (io::debug>= 2) {
                        cout << "Reference " << s.id << "\t" << a << endl;
                    }

                    if (a.compare("?")==0 || a.compare (".")==0 || a.compare("N")==0 || a.compare("-")==0) {
                        cout << "Filtering " << s.id << " : no reference base" << endl;
                        continue;
                    } else { 
                        if (useonlyconfidentsites && islower(a[0])) {
                            cout << "Filtering " << s.id << " : not a  confident site" << endl;
                            continue;
                        } else
                            boost::to_upper(a);

                        if (s.ref == 'X' && s.var=='X') { 
                            if (s.var==a[0]){
                                s.var = s.ref;
                                s.ref = a[0];
                                if (isgenotype) {
                                    for (int l = 0 ; l < nind; l++)
                                        gtype[l] = 2 - gtype[l];
                                } else {
                                    for (int l = 0 ; l < nind; l++)
                                        gtype[l] = 1 - gtype[l];
                                }
                            } else if (s.ref==a[0]) {}
                            else {
                                // Triallelic
                                cout << "Filtering " << s.id << " : triallelic" << endl;
                                continue;
                            }
                        } else  {
                            if (s.ref=='X') {
                                if (s.var==a[0]){
                                    s.ref = a[0];
                                    s.var = 'X';
                                    if (isgenotype) { 
                                        for (int l = 0 ; l < nind; l++)
                                            gtype[l] = 2 - gtype[l];
                                    } else {
                                        for (int l = 0 ; l < nind; l++)
                                            gtype[l] = 1 - gtype[l];
                                    }
                                } else {
                                    s.ref = a[0];
                                }
                            } else if (s.var=='X'){
                                if (s.ref==a[0]) {
                                } else {
                                    s.var = s.ref;
                                    s.ref = a[0];
                                    if (isgenotype) {
                                        for (int l = 0 ; l < nind; l++)
                                            gtype[l] = 2 - gtype[l];
                                    } else { 
                                        for (int l = 0 ; l < nind; l++)
                                            gtype[l] = 1 - gtype[l];
                                    }
                                }
                            }
                        }
                    }
                } else  {
                    cout << "Filtering " << s.id << " : site not in reference" << endl;
                    continue;
                }
			} 

			if (filter) { 
				if (filtermapsnps->find (s.id) != filtermapsnps->end()){
					string v = (*filtermapsnps)[s.id];
					vector<string> toks;
					functions::tokenize (v.c_str(), toks, " ");
					string a = string(1,toks[1][0]);
					if (a.compare("?")==0)
						continue;
					else { 
						if ( s.ref=='X' && s.var=='X'){
							if (s.ref==a[0] && s.var==a[0])
								continue;
						}
					}
				} else 
					continue;
			}
			if (uniqpos && flag && s.genpos<=prevgpos)
					continue;
			prevgpos = s.genpos;

			flag = true;
			vector<int> &gtype = s.gtype;

            if (outputformat.compare ("eigenstrat")==0){ 
                if (newindivname) {
                    if (keeporder) {
                        bool flag = false;
                        for (int l = 0 ; l < nind; l++)  {
                            string id = g->indiv[l].id; 
                            if (newindmap->find(id)!=newindmap->end()) {
                                genofs << gtype[l];
                                flag = true;
                            }
                        }
                        if (flag ){
                            genofs<<endl;
                        }
                    } else {
                        bool flag = false;
                        for (int l = 0; l< newindvector->size(); l++){
                            string id = (*newindvector)[l];
                            if ((g->indivmap).find(id) != (g->indivmap).end()){
                                flag = true;
                                int l1 = g->indivmap[id];
                                genofs << gtype[l1];
                            }
                        }
                        if (flag ){
                            genofs<<endl;
                        }
                    }
                } else if (extracttrios) {
                    if (!isgenotype) { 
                        unordered_map <string, int> a;
                        unordered_map <string, int> b;
                        unordered_map <string, int> counter;
                        for (int l = 0 ; l < nind/2; l++)  {
                            string id = g->indiv[2*l].id; 
                            string id1  = id.substr(0,id.length()-1);
                            char c = tolower(id[id.length()-1]);
                             if (c=='a') {
                                a[id1] = l;
                             } else if (c=='b') {
                                b[id1] = l;
                             }
                            if (counter.find(id1)!=counter.end())
                                counter[id1]=counter[id1]+1;
                             else 
                                counter[id1] = 1;

                        }

                        bool inconsistent = false;
                        for (int l = 0 ;  l < nind/2; l++){ 
                            string id = g->indiv[2*l].id; 
                            string id1  = id.substr(0,id.length()-1);
                            char c = tolower(id[id.length()-1]);
                            int aind = a[id1];
                            int bind = b[id1];
                            if ( c=='c' || c== 'd' ) {
                                if (counter[id1]==3){
                                    int a1 = gtype[2*aind]+gtype[2*aind+1]; a1-=gtype[2*l];
                                    int b1 = gtype[2*bind]+gtype[2*bind+1]; b1-=gtype[2*l+1];
                                    if (a1<0 || b1<0 || a1>1 || b1>1)
                                        inconsistent = true;

                                }
                            } 
                        }

                        if (inconsistent) { 
                            cerr << "Filtering " << s.id << ": Inconsistent genotypes" << endl;
                            continue;
                        }

                        for (int l = 0 ;  l < nind/2; l++){ 
                            string id = g->indiv[2*l].id; 
                            string id1  = id.substr(0,id.length()-1);
                            char c = tolower(id[id.length()-1]);
                            int aind = a[id1];
                            int bind = b[id1];
                            if ( c=='c' || c== 'd' ) {
                                if (counter[id1]==3){
                                    int a1 = gtype[2*aind]+gtype[2*aind+1]; a1-=gtype[2*l];
                                    int b1 = gtype[2*bind]+gtype[2*bind+1]; b1-=gtype[2*l+1];
                                    // fill
                                    genofs << gtype[2*l] << gtype[2*l+1];
                                    genofs << a1 << b1;
                                    if (trioflag) {
                                        string id2 = id + "u";
                                        trioind.push_back (id);
                                        trioind.push_back (id);
                                        trioind.push_back (id2);
                                        trioind.push_back (id2);
                                    }
                                }
                            }
                        }
                        genofs << endl;
                        trioflag = false;
 
                    }
                } else {
                    for (int l = 0 ; l < nind; l++)  {
                        genofs << gtype[l];
                    }
                    genofs<<endl;
                }
            } else if (outputformat.compare ("ihs")==0) { 
                for (int  l = 0 ; l < nind; l++){
                    outputgeno[l][outputsnps] = 1 - gtype[l];
                }
                outputsnps++;
            } else if (outputformat.compare ("beagle") ==0 ){

            }

            if (outputformat.compare ("eigenstrat")==0)
    			snpfs << s.id << "\t" << s.chr <<"\t"<<s.genpos<<"\t"<<s.getphyspos() <<"\t"<<s.var<<"\t"<<s.ref<<endl;
            else if (outputformat.compare ("ihs")==0)
    			snpfs << s.id << "\t" << s.chr <<"\t"<<s.genpos<<"\t"<<s.getphyspos() <<"\t"<<s.var<<"\t"<<s.ref<<endl;

		}
	}
	snpfs.close();
    if (outputformat.compare ("ihs")==0) {
        for (int l  = 0 ; l < nind; l++) {
            for (int i = 0 ; i<outputsnps; i++) 
                genofs << outputgeno[l][i] << "\t";
            genofs <<endl;
        }
    }

    genofs.close();
    if (extracttrios) { 
        for (int l = 0 ;  l < trioind.size(); l++) { 
            indivfs << trioind[l] <<"\tU\tU"<<endl; 
        }
    } else {
        if (keeporder){
            for (int l = 0 ; l < nind; l++){
                ind &i = g->indiv[l];
                if (newindivname) {
                    if (newindmap->find(i.id)!=newindmap->end())
                        indivfs << i.id <<"\t" <<i.gender <<"\t"<<i.group<<endl; 
                } else {
                    indivfs << i.id <<"\t" <<i.gender <<"\t"<<i.group<<endl; 
                }
            }
        } else { 
            for (int l = 0;  l < newindvector->size(); l++){
                string id = (*newindvector)[l];
                if ((g->indivmap).find(id)!=(g->indivmap).end()){
                    int l1 = g->indivmap[id];
                    ind &i = g->indiv[l1];
                    indivfs << i.id <<"\t" <<i.gender <<"\t"<<i.group<<endl; 
                }
            }
        }
    }
	indivfs.close();

}


int main (int argc, char* argv[]) {
	convertf c = convertf(argc,argv);
	c.convert ();
}
