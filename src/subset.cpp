#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "subset.h"

ofstream outfs;

subset::subset (int argc, char *argv[]) {
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
    g = new genotype (snp1,ind1, geno1, inputformat, isgenotype, "", "" );

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

void subset::convert ()  {

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
	bool pickchr = p->get_string ("chrom",chr,false);
    bool keeporder  = true;
    p->get_boolean ("keeporder", keeporder, true, false);
	bool badsnps = p->get_string ("badsnpname",gsfile,false);
    bool goodsnpstart = p->get_string ("goodsnpstart",gstart, false);
	p->get_boolean ("uniqpos", uniqpos, false, false);

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


	if (goodsnps && badsnps) {
		cerr << "Can specify only one of goodsnpname or badsnpname" <<endl;
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


	ofstream snpfs (snpoutfile.c_str());
	ofstream indivfs (indivoutfile.c_str());
	ofstream genofs (genotypeoutfile.c_str());


	int nind = g->nind;
	int nchr = g->nchr;
    int nsnps = g->nsnp;
    vector< vector <int> > outputgeno;
    int outputsnps = 0 ;

    cout << nind << "\t" << nchr << "\t" << nsnps << endl;

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

           
			if (uniqpos && flag && s.genpos<=prevgpos)
					continue;
			prevgpos = s.genpos;

			flag = true;

            if (outputformat.compare ("eigenstrat")==0){ 
                for (int l = 0 ; l < nind; l++)  {
                    int geno  = (*g)(i, j, l);
                    genofs << geno;
                }
                genofs<<endl;
            }

            if (outputformat.compare ("eigenstrat")==0)
    			snpfs << s.id << "\t" << s.chr <<"\t"<<s.genpos<<"\t"<<s.getphyspos() <<"\t"<<s.var<<"\t"<<s.ref<<endl;
            else if (outputformat.compare ("ihs")==0)
    			snpfs << s.id << "\t" << s.chr <<"\t"<<s.genpos<<"\t"<<s.getphyspos() <<"\t"<<s.var<<"\t"<<s.ref<<endl;

		}
	}
	snpfs.close();

    genofs.close();
    if (keeporder){
        for (int l = 0 ; l < nind; l++){
            ind &i = g->indiv[l];
            indivfs << i.id <<"\t" <<i.gender <<"\t"<<i.group<<endl; 
        }
    }



    indivfs.close();

}


int main (int argc, char* argv[]) {
	subset c = subset(argc,argv);
	c.convert ();
}
