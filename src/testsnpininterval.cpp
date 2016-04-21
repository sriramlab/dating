#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "convertf.h"
#include "intervalmap.h"

class testsnpininterval  {
    public:
        testsnpininterval (int argc, char *argv[]);
        void set_intervals ( ) ;
        double statistic ( unordered_map<string,string> &testset) ; 
        pair<double,double> pvalue ( double observed ) ;

        snpmap *bgsnps;
        intervalmap *imap;
        unordered_map<string, double> freqmap;
        bool givenfreq;
        vector<vector<string> > freqbins;
        vector<int> targetfreq;
        unordered_map<string,string> targetsnps;

	    unordered_map<string,vector<int> > intspersnp ;
	    unordered_map<string,int > countspersnp ;
        int iters;

        bool usemaf ;

        int bins;

        // Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

        int seed;
		const gsl_rng_type * rng_T;
        gsl_rng * rng_r;

};



testsnpininterval::testsnpininterval (int argc, char *argv[]) {
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

    p->get_int ("seed",seed,1,false);

    gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	rng_r = gsl_rng_alloc (rng_T);
	gsl_rng_set (rng_r, seed);
	p->get_int ("debug", io::debug, 0);
	p->get_int ("iters", iters, 100);


    string backgroundsnpfile;
    p->get_string ("backgroundsnps", backgroundsnpfile, true);
    bgsnps = new snpmap (backgroundsnpfile);
    for (unordered_set<snp,boost::hash<snp> >::iterator i  = bgsnps->allsnps.begin(); i != bgsnps->allsnps.end(); i++){ 
        countspersnp[i->id] = 0;
    }
    string intervalfile;
    p->get_string ("intervals",intervalfile, true);
    imap = new intervalmap (intervalfile, rng_r, bgsnps);


    string freqfile ; 
    givenfreq =  p->get_string ("freqfile", freqfile, false);
    p->get_int ("bins", bins, 10);


	p->get_boolean ("usemaf", usemaf, false);

    
    if (givenfreq) { 
	    unordered_map<string, string> *tmpmap 
             = new unordered_map<string,string> ();
        freqbins.resize (bins+1, vector<string>());
        fileio::read_map (freqfile, tmpmap, 0, 1);
    	for (unordered_map<string,string>::iterator i = tmpmap->begin(); i != tmpmap->end(); i++) { 
            string id = i->first;
            double f =  atof ((i->second).c_str());
            if (usemaf && f > 0.5)
                f=  1- f;
		    freqmap[id] = f;
            int bin = round(bins*f);
            freqbins [bin].push_back (id);
	    }

        if (io::debug >=1 ){
            cout << "Frequency bins " << endl;
            for (int i  = 0 ; i < freqbins.size() ; i++) {
                cout <<i << "\t" << freqbins[i].size() << endl;
            }
        } else if (io::debug>=2) { 
            cout << "Frequency bins " << endl;
            for (int i  = 0 ; i < freqbins.size(); i++) {
                cout <<i << "\t" << freqbins[i].size() << endl;
                for (int j = 0 ; j < freqbins[i].size(); j++)
                    cout << freqbins[i][j] << "\t";
                cout <<endl;
            }
        }
    }

    string targetsnpfile;
    bool targetset = p->get_string ("targetsnps", targetsnpfile, true);
    fileio::read_map (targetsnpfile, &targetsnps, 0);
    if (givenfreq) { 
        targetfreq.resize(bins+1);
        for ( unordered_map<string, string> ::iterator i = targetsnps.begin(); i != targetsnps.end(); i++){
            string id = i->first;
            if (freqmap.find(id) ==freqmap.end()) {
                cerr  << "Could not find " << id << " in the background file" <<endl;
                exit(1);
            }
            double f = freqmap[id];
            int bin = round(bins * f);
            targetfreq[bin]++;
        }   

        if (io::debug >=1 ){
            cout << "Target frequency bins " << endl;
            for (int i  = 0 ; i < targetfreq.size() ; i++) {
                cout <<i << "\t" << targetfreq[i] << endl;
            }
        }
    }
    set_intervals ();
  
    string iefile;
    bool intervalenrichment =  p->get_string ("enrichmentstatistic",iefile, false);
    if (intervalenrichment) { 
  	    unordered_map<string, string> *tmpmap 
             = new unordered_map<string,string> ();
        fileio::read_map (iefile, tmpmap, 0, -1);
        unordered_map <string,double> emap;
        unordered_map <string,int> smap;
    	for (unordered_map<string,string>::iterator i = tmpmap->begin(); i != tmpmap->end(); i++) { 
            string id = i->first;
            string v = i->second;
            vector<string> toks;

            functions::tokenize (v.c_str(), toks, "\t ");

            double val1 = atof (toks[1].c_str());
            emap[id] =val1;
            int val2 = atoi(toks[2].c_str());
            smap[id] = val2;
        }
	    unordered_map<string, vector<interval> > &intervals = imap->intervals;
    	unordered_map <string, vector<interval> >::iterator it;

        ofstream efs ("enrichment.txt");
        for (it = intervals.begin(); it != intervals.end (); it++) { 
            vector<interval> & tmpints = it->second;
            string chr = it->first;
		    vector<snp> &snps = (bgsnps->snps)[chr];
    		for (int j = 0 ; j < tmpints.size(); j++){ 
	    		interval &i = tmpints[j];

    			vector<int> &snpind = (imap->snpsperint)[i];
                double avgstat1 = 0 ;
                double avgstat2 = 0 ;
                int denom =  0;
			    for (int k = 0 ; k < snpind.size(); k++){ 
				    snp &snp = snps[snpind[k]];
                    avgstat1 += emap[snp.id];
                    avgstat2 += smap[snp.id];
                    denom++;
                }
                if (denom > 0 ) {
                    avgstat1/=denom;
                    avgstat2/=denom;
                }
                efs << i.to_string() << "\t" << avgstat1 << "\t" << avgstat2 << "\t" << denom << endl;
            }
        }
        int nchr = bgsnps->nchr;

        ofstream nefs ("non-enrichment.txt");
        for (int i =  0; i  <nchr; i ++) {
            string chr = (bgsnps->chrs)[i];
            vector<snp> &snps = (bgsnps->snps)[chr];
            for (int j = 0 ; j < snps.size(); j++){
                snp &s = snps[j];
               if (intspersnp[s.id].size() == 0  ) {
                    double avgstat1 = emap[s.id];
                    double avgstat2 = smap[s.id];
                    nefs << avgstat1 << "\t" << avgstat2 << endl;
               }
            }
        }
        nefs.close ();

        efs.close ();
        exit(0);
    }


    if (io::debug>=1 ) {
        cout << "Target snps " << endl;
        for ( unordered_map<string, string> ::iterator i = targetsnps.begin(); i != targetsnps.end(); i++){
            string id =  i->first;
            vector<int>& tmp = intspersnp[id];
            cout << id << "\t" << tmp.size() <<"\t" << countspersnp[id] << endl;
        } 
    }
} 

double testsnpininterval::statistic ( unordered_map<string,string> &testset) { 
    int total = 0 ;
    for (unordered_map<string,string>::iterator i =  testset.begin(); i != testset.end(); i++) { 
            string id =  i->first;
            total += (countspersnp[id] > 0);
    }
    return (total);
}



pair<double,double> testsnpininterval::pvalue ( double observed ) {
    double pval1 = 0 ;
    double pval2 = 0 ;
    vector<int> tmpfreq;
    tmpfreq.resize(bins+1);
    for (int i  = 0  ; i < iters; i ++){ 
        cout << "iter "  << i << endl;
        if (givenfreq) { 
            int stat =  0;
            for (int j = 0 ; j < tmpfreq.size(); j++)
                tmpfreq[j] = 0 ;

            for (int j = 0 ; j < freqbins.size(); j ++){
                int k = targetfreq[j];
                int n = freqbins[j].size();
                if (n==0)
                    continue;
                int* base = new int[n];
                vector<string> & tmp = freqbins[j];
                for (int l = 0 ; l  < n; l++)
                    base[l] = l;
                gsl_ran_shuffle (rng_r, base, n, sizeof(int));

               if (io::debug>=1) { 
                   cout << "k = " << k << endl;
                   if (io::debug >= 2 ) {
                       cout << "Permuted snp ids " <<endl;
                       for (int l  = 0 ; l < tmp.size(); l++)
                           cout << tmp[l] << "\t";
                       cout << endl;
                       for (int l  = 0 ; l < tmp.size(); l++)
                           cout << tmp[base[l]] << "\t";
                       cout << endl;
                   }
               }


                int l = 0  ;
                int index =  0 ;
                vector<string> ids;
                while (index < n) {
                    string id = tmp[base[index]];
                    if (l >= k)
                        break;
                    if (targetsnps.find (id)==targetsnps.end()){
                        ids.push_back(id);
                        stat += (countspersnp[id] > 0 );           
                        l++;
                    }
                    index++;
                }

                if (io::debug>=1) {

                    cout << "Chosen ids"<<endl;
                    for (int i  = 0 ; i < ids.size();i++)
                        cout << ids[i] << "\t";
                    cout << endl;
                    tmpfreq[j] += ids.size();
                }

                if (index==n-1){
                    cerr << "Could not match frequencies in bin " << j << endl;
                    exit (1);
                }
                delete[] base;
            }
            if (io::debug >= 1){
                cout << "Tmp frequency bins " << endl;
                for (int i  = 0 ; i < tmpfreq.size() ; i++) {
                    cout <<i << "\t" << tmpfreq[i] << endl;
                }
            }

            pval1 += (stat >= observed);
            pval2 += (stat <= observed);
            cout << "Iteration: " << i << "\t" << stat << "\t" << observed << endl;
        }
         
    }
    if (iters>0)
        pval1 /= iters;
    if (iters>0)
        pval2 /= iters;
    return pair<double,double>(pval1,pval2);
}

void testsnpininterval::set_intervals ( ) {
	unordered_map<interval, vector<int>, boost::hash<interval> > & snpsperint = imap->snpsperint;
	unordered_map<string, vector<interval> > &intervals = imap->intervals;
	unordered_map<string, vector<snp> >& snps = bgsnps->snps;	
	double eps = 1e-9;

	io::println ("In set_intervals", 2);
	unordered_map <string, vector<interval> >::iterator it;
	for (it = intervals.begin(); it != intervals.end (); it++) { 
		vector<interval> & tmpints = it->second;
		string chr = it->first;
		if (snps.find(chr)==snps.end())
			continue;
		vector<snp> & tmpsnps = snps[chr];
		int j = 0;
		for (int i = 0 ; i < tmpints.size(); i++) {
			double s =  tmpints[i].start;
			double e = tmpints[i].end;
			io::println ("Examining interval " + tmpints[i].to_string(),2);
			while ( j < tmpsnps.size() && tmpsnps[j].physpos < s)
				j++;

			if (j >= tmpsnps.size())
				break;
			int k = j;
			while ( k < tmpsnps.size() && tmpsnps[k].physpos <= e) {
				if (!tmpsnps[k].ignore) {
                    string id = tmpsnps[k].id;
					intspersnp[id].push_back(i);
                    countspersnp[id]++;
					snpsperint[tmpints[i]].push_back(k);

					snp& ts = tmpsnps[k];
					ts.insideinterval =  true;
					if ( abs(ts.physpos - s) <eps )
						ts.start = true;
					if ( abs(ts.physpos - e) <eps )
						ts.end = true;
					io::println ("Setting  SNP "  + ts.to_string(),2);
				}
				k++;
			}
		}
	}

}

int main (int argc, char* argv[]) {

    testsnpininterval test= testsnpininterval (argc, argv);
    double observed = test.statistic (test.targetsnps);
    cout << "Observed =  " << observed << endl;
    pair<double,double> pvpair = test.pvalue (observed);
    cout << "Pval = " << pvpair.first << "\t" << pvpair.second << endl;
}
