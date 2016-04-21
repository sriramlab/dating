#include "data.h"
#include "fileio.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "computed.h"
#include<boost/filesystem/operations.hpp>
namespace bf = boost::filesystem;		 //create an alias

ofstream outfs;
ofstream pairfs;

computed::computed (int argc, char *argv[]) {
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
	out = p->get_string ("output",outfile,true);
	if (out) {
	       outfs.open (outfile.c_str());	
	}

	pairs = p->get_string ("pairfile", pairfile, false);
	toflip  =p->get_string ("flipsnpname", flipfile, false);
	if (toflip) {
		fileio::read_map (flipfile, &flipmap, 0);
	}
	p->get_double ("maxdis",maxdis,0.01);
	p->get_double ("binsize",binsize,0.01);
	p->get_int ("debug", io::debug, 0);
    p->get_boolean ("sample", sample, false, false);
	goodsnps = p->get_string ("goodsnpname",gsfile,false);
	goodsnpsmap = NULL;
	if (goodsnps ){
		goodsnpsmap = new unordered_map <string, string > ();
		fileio::read_map (gsfile, goodsnpsmap,0);
		if (io::debug >= 2)  {
			typedef unordered_map<string, string>::iterator iter;
			for ( iter i  = goodsnpsmap->begin();  i != goodsnpsmap->end(); i++){
				cout << i->first << "\t" << i->second << endl;
			}
		}
	}
	badsnps = p->get_string ("badsnpname",gsfile,false);
	badsnpsmap = NULL;
	p->get_int ("ascertain", ascertain, 0);
    p->get_boolean ("isgenotype", isgenotype, true);
	if (badsnps){
		badsnpsmap = new unordered_map <string, string > ();
		fileio::read_map (gsfile, badsnpsmap,0);
		if (io::debug >= 2)  {
			typedef unordered_map<string, string>::iterator iter;
			for ( iter i  = badsnpsmap->begin();  i != badsnpsmap->end(); i++){
				cout << i->first << "\t" << i->second << endl;
			}
		}
	}

    filtersnps = "none";
    givensnpfilter = p->get_string ("filtersnps", filtersnps, false);

	if (pairs)  {
		string weightfile;
		p->get_string ("weightname",weightfile,true);
		read_weights (weightfile);
	} else { 
		g = new genotype (p);
		outputpair = p->get_string ("outputpair",outputpairfile,false);
		dojack = p->get_string ("jackknife", jackpath, false);
        jackblock = p->get_double ("jackblock", jackblock, -1, false);
        jackblock = 1e6 * jackblock;
        p->get_boolean ("jackgenetic", jackgenetic, false, false);
        p->get_int ("jackminsnps", jackminsnps, 100, false);
		if (outputpair) { 
			pairfs.open (outputpairfile.c_str());
		}
	}
}

int computed::positiontojbin (int chr, int index) { 
    if (jackblock <= 0){
        cerr << "Invalid jackknife block size" << endl;
        exit(1);
    }
    vector<snp> &tmpsnps =  g->snps[g->chrs[chr]];
    int pos = tmpsnps[index].physpos;
    int i = (pos - g->minppos[g->chrs[chr]])/jackblock;
    i = i + cumnumjbins[chr];
    return i;    
}

void computed::compute () {
	nbins =  mathfn::round (maxdis/binsize );
	vector<double> ld1 (nbins);
	vector<int> lddenom1 (nbins);
	int nsnp1 = 0;
	int nind = g->nind;
	int nchr = g->nchr;
    int samplecount = 0 ;
    
    int njbins = nchr;
	vector<vector<double> > ldj  (nchr, vector<double>(nbins));
	vector<vector<int> > lddenomj(nchr, vector<int>(nbins));
    cumnumjbins.resize (nchr, 0);
    if (jackblock > 0 ){
        njbins =  0;
        for (int i = 0 ; i < nchr ; i++){
            int tmp = ceil( (g->maxppos[g->chrs[i]] - g->minppos[g->chrs[i]])/jackblock) ;
            cumnumjbins[i] = njbins;
            njbins += tmp;
        }
        ldj.resize (njbins, vector<double> (nbins));
        lddenomj.resize (njbins, vector<int> (nbins));
    } else {
        njbins = nchr;
    } 
    vector<int> g1 (g->nind, 0);
    vector<int> g2 (g->nind, 0);
	for (int i = 0 ; i < nchr; i++){  
		vector<snp> &tmpsnps =  g->snps[g->chrs[i]];
		io::println ("chr = " + tostring (i), 0);
		vector<double> &tmpld = ldj[i];
		vector<int> &tmplddenom = lddenomj[i];

		for (int j  = 0 ; j < tmpsnps.size();j++) {
			for (int k  = j+1; k < tmpsnps.size (); k++) {
				snp &s1 = tmpsnps[j];
				snp &s2 = tmpsnps[k];
                g->get_geno (i,j,g1);
                g->get_geno (i,k,g2);
//				vector<int>& g1 = s1.gtype;
//				vector<int>& g2 = s2.gtype;
                if (filtersnps.compare ("transitions")==0 || filtersnps.compare("transition")==0) {
                    if (s1.istransition() && s2.istransition()) {}
                    else continue;
                } else if (filtersnps.compare("transversions")==0 || filtersnps.compare("transversion")==0) {
                    if (s1.istransversion() && s2.istransversion()) {}
                    else continue;
                }

                if (goodsnps){
                    if (goodsnpsmap->find(s1.id) == goodsnpsmap->end() || goodsnpsmap->find(s2.id) == goodsnpsmap->end())
                        continue;
                } 
                if (badsnps)  {
                    if (badsnpsmap->find(s1.id) != badsnpsmap->end() || badsnpsmap->find(s2.id) != badsnpsmap->end())
                        continue;
                }

				double dis = s2.genpos - s1.genpos;
				double pdis = s2.physpos - s1.physpos;

				if (s2.genpos  - s1.genpos > maxdis)
					continue;
				int bini  = (int)(dis/binsize);

				double c[3][3];
				int denom = 0 ;
				for (int i1 =  0 ; i1 < 3; i1++)
					for (int j1= 0 ; j1 < 3; j1++)
						c[i1][j1] = 0;

				for (int l = 0 ; l < nind; l++){
					io::println ( "g1 =  " + tostring (g1[l] ) + "\t g2 = " + tostring(g2[l]),2);
					if ( g1[l]==9 || g2[l]==9)
						continue;
					c[g1[l]][g2[l]]++;
					denom++;
				}

				double d1 = 0 ;
				double p1 = 0 ; 
				double q1 = 0;
				for (int i1 = 0 ; i1 < 3; i1 ++){ 
					for (int j1  = 0; j1 < 3; j1++){ 
						d1 += i1*j1*c[i1][j1];
						p1 += i1*c[i1][j1];
						q1 += j1*c[i1][j1];
					}
				}


				d1/=denom;
				d1-= p1*q1/(denom*denom);
				double r2 = 0;
                
                int factor = 1;
                if (isgenotype==1)
                    factor = 2;

				p1 /= (factor*denom);
				q1 /= (factor*denom);
				if (p1 > 0 && p1 < 1 && q1 > 0 && q1 < 1)
					r2 = (d1*d1)/(factor*factor*p1*(1-p1)*q1*(1-q1));

				if (outputpair) { 
					pairfs <<  s1.chr << "\t" <<  s1.getphyspos() << "\t" << s2.chr << "\t" << s2.getphyspos() << "\t" << dis << "\t" << pdis << "\t" << d1 << "\t" << p1 << "\t" << q1 << "\t" << r2 << endl;
				}

				ld1[bini] += d1;
				lddenom1[bini] ++;
				if (dojack) {
                    if (jackblock >0 ) {
                        int jj = positiontojbin (i, j);
                        int jk = positiontojbin (i, k);

                        for (int ji = jj ; ji <= jk; ji ++ ){ 
                            vector<double> &tmpld = ldj[ji];
                            vector<int> &tmplddenom = lddenomj[ji];
                            tmpld [bini] += d1;
                            tmplddenom [bini] ++;

                        }
                    } else {
					    tmpld[bini] += d1;
					    tmplddenom[bini] ++;
                    }
				}
				io::println ("d1 = "  + tostring(d1),2);

                if (sample) {
                    cout.setf(ios::fixed,ios::floatfield); 
	                cout.precision(7);
					cout <<  s1.chr << "\t" <<  s1.getphyspos() << "\t" << s2.chr << "\t" << s2.getphyspos() << "\t" <<  d1 << "\t" << p1 << "\t" << q1 <<  endl;
                    samplecount ++;
                    if (samplecount  > 10000)
                        exit(1);
                }


			}
		}
	}
	for (int i = 0 ; i < nbins; i++) {
		if (lddenom1[i]>0){
			ld1[i]/=lddenom1[i];
		}
//		double j = (i+1.0)/nbins;
		double j = (i+1.0)*binsize*100;
		if (out) { 
			outfs << j <<"\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << lddenom1[i] ;
			outfs << endl;
		} else { 
			cout << j <<"\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << lddenom1[i] ;
			cout << endl;
		}	
	}

	if (dojack) {
		bf::create_directory (bf::path (jackpath));
		bf::path p (jackpath);
		p /= "snps.txt";
		ofstream weightfs (p.string().c_str());
		for (int  i = 0, j = 0 ; i < njbins; i++){ 
			vector<double> &tmpld = ldj[i];
			vector<int> &tmplddenom = lddenomj[i];
			int snps =  0;
			for (int k = 0 ; k < nbins; k++) {
				snps += tmplddenom[k];
            }
            if (snps  < jackminsnps) 
                continue;
            j++;

			bf::path p (jackpath);
			string fname ;
			if (out) {
				fname = outfile;
			       	fname = fname + ".";
			        fname = fname + tostring(j+1);
			} else  {
				fname = "out.";
				fname = fname + tostring(j+1);
			}
			p /= fname;

			ofstream tmpofs (p.string().c_str());

			for (int k = 0 ; k < nbins; k++) {
				tmpld[k] = ld1[k] * lddenom1[k] - tmpld[k];
				tmplddenom[k] = lddenom1[k] - tmplddenom[k];
				if (tmplddenom[k]>0){
					tmpld[k]/=tmplddenom[k];
				}
				double j = (k+1.0)/nbins;
				tmpofs << j <<"\t" << tmpld[k] << "\t" << tmpld[k] << "\t" << tmpld[k] << "\t" << tmpld[k] << "\t" << tmplddenom[k] ;
				tmpofs << endl;

			}
			weightfs << snps << endl;
			tmpofs.close ();
		}
		weightfs.close ();
	}

	outfs.close ();
}

void computed::read_weights (string filename) { 
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


void computed::computepairs () { 
	nbins =  mathfn::round (maxdis/binsize );
	vector<double> ld1 (nbins);
	vector<int> lddenom1 (nbins);
	int nsnp1 = 0;
	ifstream inp (pairfile.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< pairfile <<endl;
		exit(1);
	}
	io::println ("In computepairs",1);
	string line;
	int linenum = 0;
    int samplecount = 0;
	while ( std::getline (inp, line) ) {
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);

		string s1;
		string s2;
		string pos1;
		string pos2;
		string id1;
		string id2;
		double pdis;
		double gdis;
		double p[4];
		double d;
		ss >> s1; 
		ss >> pos1;
		ss >> s2;
		ss >> pos2;
		ss >> gdis;
		ss >> pdis;

		// 0=ref-ref
		// 1=ref-alt
		// 2=alt-ref
		// 3=alt-alt
		for (int i  = 0 ; i < 4; i++)
			ss >> p[i];


		id1 = "chr" + s1 + ":" + pos1;
		id2 = "chr" + s2 + ":" + pos2;
		if (goodsnps){
			if (goodsnpsmap->find(id1) == goodsnpsmap->end() || goodsnpsmap->find(id2) == goodsnpsmap->end())
				continue;
		} 
		if (badsnps)  {
			if (badsnpsmap->find(id1) != badsnpsmap->end() || badsnpsmap->find(id2) != badsnpsmap->end())
				continue;
		}

		io::println (tostring (id1) + "\t" + tostring(id2) + "\t"+ tostring(gdis) + "\t" + tostring(pdis) + "\t" +  tostring (p[0]) + "," + tostring(p[1]) + "," + tostring(p[2]) + "," + tostring(p[3]),1);
		
		if ( weightmap.find (id1) == weightmap.end () || weightmap.find(id2) == weightmap.end())
			continue;
		double w1  = weightmap[id1];
		double w2 = weightmap[id2];

		if (flipmap.find (id1) != flipmap.end())
			w1 = 1-w1;
		if (flipmap.find (id2) != flipmap.end())
			w2 = 1-w2;

		if (ascertain ==1 ) { 
			if (w1<=0.1 && w2<=0.1) {}
			else continue;

		}

		if (flipmap.find(id1)!=flipmap.end()){
			double tmp = p[2];
			p[2] = p[0];
			p[0] = tmp;

			tmp = p[1];
			p[1] = p[3];
			p[3] = tmp;
		}
		if (flipmap.find(id2)!=flipmap.end()){
			double tmp = p[0];
			p[0] = p[1];
			p[1] = tmp;

			tmp = p[2];
			p[2] = p[3];
			p[3] = tmp;
		}

		io::println (tostring (id1) + "\t" + tostring(id2) + "\t"+ tostring(gdis) + "\t" + tostring(pdis) + "\t" +  tostring (p[0]) + "," + tostring(p[1]) + "," + tostring(p[2]) + "," + tostring(p[3]),1);


		double dis = gdis;
		if (dis > maxdis)
			continue;
		int bini  = (int)(dis/binsize);
		double d1  = p[0] * p[3] - p[1] * p[2];

        if (sample) {
            double a = p[2]+p[3];
            double b = p[1]+p[3];
            cout.setf(ios::fixed,ios::floatfield); 
	        cout.precision(7);
            cout << s1 << "\t" << pos1 << "\t" << s2 << "\t" << pos2 << "\t" << d1 <<"\t" << a <<"\t"<<b<<endl;
            samplecount++;
            if (samplecount > 10000)
                exit(1);
        }

		ld1[bini] += d1;
		lddenom1[bini] ++;
		io::println ("d1 = "  + tostring(d1),1);
	}	
	inp.close ();

	for (int i = 0 ; i < nbins; i++) {
		if (lddenom1[i]>0){
			ld1[i]/=lddenom1[i];
		}
		double j = (i+1.0)/nbins;
		if (out) { 
			outfs << j <<"\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << lddenom1[i] ;
			outfs << endl;
		} else { 
			cout << j <<"\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << ld1[i] << "\t" << lddenom1[i] ;
			cout << endl;
		}	
	}

	outfs.close ();

}

int main (int argc, char *argv[]) {
	computed cd (argc, argv);
	if (cd.pairs) 
		cd.computepairs ();
	else 
		cd.compute ();
}
