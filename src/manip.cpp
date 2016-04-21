#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "manip.h"

ofstream outfs;

manip::manip (int argc, char *argv[]) {
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
	p = new data (pfile, d->configMap);
	p->print_parameters();
	
	string geno1, snp1, ind1;
    inputformat = "eigenstrat";
    p->get_string ("inputformat", inputformat, false);
	p->get_string ("genotypename",geno1,true);
	p->get_string ("snpname",snp1,true);
	p->get_string ("indivname",ind1,true);
    outputformat = "eigenstrat";
    p->get_string ("outputformat", outputformat, false);


    if (inputformat.compare ("prob")==0 || outputformat.compare("prob")==0) {
        if (inputformat.compare (outputformat)!=0){
            cerr << "Input and output format not compatible" << endl;
            exit(1);
        }
    }

	p->get_boolean ("isgenotype", isgenotype, true);
    p->get_int ("debug" , io::debug, 0 );
    string goodsnpfile =  "";
    string goodindfile = "";
    bool pack = true;
	bool goodsnps = p->get_string ("goodsnpname",goodsnpfile);
    bool goodind = p->get_string ("goodindname", goodindfile);
    p->get_boolean ("pack", pack, true);
    bool hasmissing = true;
    p->get_boolean ("hasmissing", hasmissing, true);


    if (io::debug>=2) 
        cout << "Reading genotypes from " <<geno1<<endl;

//    g = new genotype (snp1,ind1, geno1, inputformat, isgenotype, goodindfile, goodsnpfile, hasmissing, pack);
    g = new genotype (p);    

    if (io::debug>=2) 
        cout << "Read genotypes from " <<geno1<<endl;
}

void manip::run ()  {

    p->get_string("snpoutname|snpoutfilename",snpoutfile,true);
	p->get_string("indivoutname|indoutfilename",indivoutfile,true);
	p->get_string("genotypeoutname|genooutfilename",genotypeoutfile,true);

    bool nomissing;
    bool onlypolymorphic;
    bool noallmissing;
	string newindivnamefile;
	string chr;
    string gstart;
    string gend;
    string gstartchr;
    double gstartpos;
    string gendchr;
    double gendpos;


    p->get_boolean ("nomissing", nomissing, false);
    p->get_boolean ("noallmissing", noallmissing, false);
    p->get_boolean ("onlypolymorphic", onlypolymorphic, false);

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

    // Flip
    bool flip = false;
    p->get_boolean ("flip", flip, false);


    // Polarize relative to individual with id polarizeid
    string polarizeid ; 
    int polarizeind = -1;
    bool polarize = p->get_string ("polarize", polarizeid);
    bool separatefile = false;
    bool useonlyconfidentsites;

	unordered_map<string, string > *pmapsnps = NULL;
    if (polarize){
        for (int l = 0 ; l < g->nind; l++){
            ind &i = g->indiv[l];
            if (polarizeid.compare(i.id)==0)   {
                polarizeind = l;
                break;
            }
        }
        if (polarizeind < 0) {
            cerr << "Id " << polarizeid << " not found" << endl;
            exit(1);
        }
    } else {
	    string polarizefile = "";
	    polarize = p->get_string ("polarizesnpname", polarizefile, false);
        p->get_boolean("useonlyconfidentsites", useonlyconfidentsites, true, false); 
        if (polarize) {
            separatefile = true;
            pmapsnps = new unordered_map <string, string> ();
            fileio::read_map (polarizefile,  pmapsnps, 0 , 1);
            if (io::debug >= 2)  {
                typedef unordered_map<string, string>::iterator iter;
                for ( iter i  = pmapsnps->begin();  i != pmapsnps->end(); i++){
                    cout << i->first << "\t" << i->second << endl;
                }
            }
        }

    }

    //
    // Select a chromosome
    bool pickchr = p->get_string ("chrom",chr,false);

    // New ids
	bool newindivname = p->get_string ("newindivname", newindivnamefile, false);
    bool keeporder  = true;
    p->get_boolean ("keeporder", keeporder, true, false);
	
    // Only retain SNPs with non-zero genetic distance to previous SNP
    p->get_boolean ("uniqpos", uniqpos, false, false);


	snpmap *smap;
	unordered_map<string, string > *newindmap = NULL;
    vector<string> *newindvector = NULL;
	

	
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
    // snps are kept if snpfilter is true
    vector<bool> *snpfilter = new vector<bool> (g->nsnp, true);
    int outputsnps = 0 ;

    if (outputformat.compare ("prob")==0) { 
        genofs.setf(ios::fixed,ios::floatfield); 
        genofs.precision(4);
    }

    
   	for (int i = 0, snpindex = 0 ; i < nchr; i++, snpindex++){
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];
        cout << "chr :"  << i << "\t" << tmpsnps.size()<< endl;
		if (pickchr && chr!=g->chrs[i])
			continue;
		double prevgpos ;
		bool flag = false;
        vector <int> gtype (g->nind, 0);
        vector <double> prob (g->nind, 0);
		for (int j  = 0 ; j < tmpsnps.size();j++) {
			snp &s = tmpsnps[j];
            if (inputformat.compare ("prob")==0) { 
                g->get_prob (i, j, prob);
            } else { 
                g->get_geno (i, j, gtype);
            }
            if (io::debug >= 2 ) {
                cout  << j << "\t" << s.id << endl;
            }
			
            if (goodsnpstart || goodsnpend) { 
                if (goodsnpstart) {
                    if (s.chr!=gstartchr || s.physpos < gstartpos)
                        (*snpfilter)[snpindex] = false;
                        continue;
                }
                if (goodsnpend) { 
                    if (s.chr!=gendchr || s.physpos > gendpos)
                        (*snpfilter)[snpindex] = false;
                        continue;
                }
            }

            if ((outputformat.compare("packedeigenstrat")==0 || inputformat.compare ("prob")) && (polarize||flip||nomissing||onlypolymorphic||noallmissing)){
                // Hard to do these operations on some file formats
                cerr << "Incompatible operations on entry of type prob" << endl;
                cerr << "Must disable polarize, flip, nomissing, onlypolymorphic, noallmissing" << endl;
                exit(1);
            }
            if (polarize) {
                if (!separatefile) { 
                    int tmp = gtype[polarizeind];
                    if (tmp==9) {
                        cout << "Filtering " << s.id << " : no reference base" << endl;
                        (*snpfilter)[snpindex] = false;
                        continue;
                    }
                    tmp = (tmp==2)?1:tmp;  // Handle diploids
                    if (tmp > 0) {
                        char a = s.var;
                        s.var = s.ref;
                        s.ref = a;
                        if (isgenotype) {
                            for (int l = 0 ; l < nind; l++) 
                                gtype[l] = gtype[l]!=9?(2 - gtype[l]):gtype[l];
                        } else {
                            for (int l = 0 ; l < nind; l++)
                                gtype[l] = gtype[l]!=9?(1 - gtype[l]):gtype[l];
                        }
                    }
                    g->set_geno (i, j, gtype);
                } else { 
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

                            if (s.ref != 'X' && s.var!='X') { 
                                if (s.var==a[0]){
                                    s.var = s.ref;
                                    s.ref = a[0];
                                    if (isgenotype) {
                                        for (int l = 0 ; l < nind; l++)
                                            gtype[l] = gtype[l]!=9 ? (2 - gtype[l]):gtype[l];
                                    } else {
                                        for (int l = 0 ; l < nind; l++)
                                            gtype[l] = gtype[l]!=9 ?(1 - gtype[l]):gtype[l];
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
                                                gtype[l] = gtype[l]!=9 ? (2 - gtype[l]):gtype[l];
                                        } else {
                                            for (int l = 0 ; l < nind; l++)
                                                gtype[l] = gtype[l]!=9 ?(1 - gtype[l]):gtype[l];
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
                                                gtype[l] = gtype[l]!=9 ? (2 - gtype[l]):gtype[l];
                                        } else { 
                                            for (int l = 0 ; l < nind; l++)
                                                gtype[l] = gtype[l]!=9 ?(1 - gtype[l]):gtype[l];
                                        }
                                    }
                                }
                            }
                        }
                    } else  {
                        cout << "Filtering " << s.id << " : site not in reference" << endl;
                        continue;
                    }
                    g->set_geno (i, j, gtype);
                }

            } 
            if (flip) { 
                for (int l = 0 ; l < nind; l++)   {
                    if (gtype[l]!=9) { 
                        if (isgenotype ) 
                            gtype[l] = 2 - gtype[l];
                        else 
                            gtype[l] = 1 - gtype[l];
                    }
                }
                g->set_geno (i, j, gtype);
            }
            if (nomissing || onlypolymorphic || noallmissing)  {
                int vc, rc ,mc;
                g->get_freq ( gtype, vc, rc, mc);
                if (nomissing && mc > 0 ) {
                    cout << "Filtering " << s.id << " : has missing data" << endl;
                    (*snpfilter)[snpindex] = false;
                    continue;
                }
                if (onlypolymorphic && (vc == 0 || rc == 0) ) {
                    cout << "Filtering " << s.id << " : monomorphic " << endl;
                    (*snpfilter)[snpindex] = false;
                    continue;
                }
                if (noallmissing && vc ==0 && rc == 0 ) {
                    cout << "Filtering " << s.id << " : has all  missing data" << endl;
                    (*snpfilter)[snpindex] = false;
                    continue;
                }
            }

			if (uniqpos && flag && s.genpos<=prevgpos)
					continue;
			prevgpos = s.genpos;

			flag = true;

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
                } else {
                    for (int l = 0 ; l < nind; l++)  {
                        genofs << gtype[l];
                    }
                    genofs<<endl;
                }
            } else if (outputformat.compare ("prob")==0){
                if (newindivname) {
                    if (keeporder) {
                        bool flag = false;
                        for (int l = 0 ; l < nind; l++)  {
                            string id = g->indiv[l].id; 
                            if (newindmap->find(id)!=newindmap->end()) {
                                genofs << prob[l] << "\t";
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
                                genofs << prob[l1] << "\t";
                            }
                        }
                        if (flag ){
                            genofs<<endl;
                        }
                    }
                } else {
                    for (int l = 0 ; l < nind; l++)  {
                        genofs << prob[l] << "\t";
                    }
                    genofs<<endl;
                }
            } else if (outputformat.compare ("packedeigenstrat")==0){
            } else if (outputformat.compare ("oxford")==0)  {
                genofs << s.chr <<" " << s.id << " " << s.getphyspos () << " " << s.var << " " << s.ref ;
                int c[3]; 
                for (int  l = 0 ; l < nind ; l ++ ){  
                    c[0] = c[1] = c[2] = 0; 
                    if (gtype[l] != 9 ) { 
                        c[gtype[l]] = 1;
                    }
                    genofs << " " << c[2] << " " << c[1] << " " << c[0];
                }
                genofs << endl;
            }

            if (outputformat.compare ("eigenstrat")==0 || outputformat.compare ("packedeigenstrat")==0 || outputformat.compare ("prob")==0)
    			snpfs << s.id << "\t" << s.chr <<"\t"<<s.getgeneticpos ()<<"\t"<<s.getphyspos() <<"\t"<<s.var<<"\t"<<s.ref<<endl;
            else if (outputformat.compare ("ihs")==0)
    			snpfs << s.id << "\t" << s.chr <<"\t"<<s.getgeneticpos ()<<"\t"<<s.getphyspos() <<"\t"<<s.var<<"\t"<<s.ref<<endl;

		}
	}
    if (outputformat.compare ("packedeigenstrat")==0) { 
        g->write_genotypes (genofs);
    }
	snpfs.close();
    genofs.close();
    if (outputformat.compare ("eigenstrat")==0 || outputformat.compare ("packedeigenstrat")==0 || outputformat.compare ("prob")==0) {
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
    } else if (outputformat.compare ("oxford")==0) { 
        indivfs << "ID_1 ID_2 missing" << endl;
        indivfs << "0 0 0"<<endl;
        if (keeporder){
            for (int l = 0 ; l < nind; l++){
                ind &i = g->indiv[l];
                if (newindivname) {
                    if (newindmap->find(i.id)!=newindmap->end())
                        indivfs << i.id <<" " <<i.id<<"\t"<<"0"<<endl; 
                } else {
                    indivfs << i.id <<" " <<i.id <<" "<<"0"<<endl; 
                }
            }
        } else { 
            for (int l = 0;  l < newindvector->size(); l++){
                string id = (*newindvector)[l];
                if ((g->indivmap).find(id)!=(g->indivmap).end()){
                    int l1 = g->indivmap[id];
                    ind &i = g->indiv[l1];
                    indivfs << i.id <<" " <<i.id <<" "<<"0"<<endl; 
                }
            }
        }

    }

	indivfs.close();

}


int main (int argc, char* argv[]) {
	manip t = manip(argc,argv);
	t.run ();
}
