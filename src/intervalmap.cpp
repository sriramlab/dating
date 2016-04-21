#include "intervalmap.h" 
#include "snpmap.h"
#include "stringfn.h"
#include "vectorfn.h"
#include "functions.h"
#include "io.h"
#include "gzstream.h"

intervalmap::intervalmap (string filename, string chrom ) {
    vector<string> tok;
    functions::tokenize(filename.c_str(), tok,".");

    this->chrom = chrom;
    if (chrom.compare ("")!=0)
        givenchrom =  true;
    else 
        givenchrom = false;

    if (tok[tok.size()-1].compare ("gz")==0 || tok[tok.size()-1].compare("bgz")==0){
        igzstream ig (filename.c_str());
        getstats(ig);
        ig.close();
        igzstream ig1 (filename.c_str());
        read_stream (ig1, r, smap);
        ig1.close ();
    } else {
        ifstream inp (filename.c_str());
        if (!inp.is_open()){
            cerr << "Error reading file "<< filename <<endl;
            exit(1);
        }
        getstats(inp);
        inp.close ();
        ifstream inp1 (filename.c_str());
        read_stream (inp1, r, smap);
        inp1.close();
    }
}

intervalmap::intervalmap (string filename, gsl_rng *r, snpmap *smap ) {
    vector<string> tok;
    functions::tokenize(filename.c_str(), tok,".");
    givenchrom = false;

    if (tok[tok.size()-1].compare ("gz")==0 || tok[tok.size()-1].compare("bgz")==0){
        igzstream ig (filename.c_str());
        getstats(ig);
        ig.close();
        igzstream ig1 (filename.c_str());
        read_stream (ig1, r, smap);
        ig1.close ();
    } else {
        ifstream inp (filename.c_str());
        if (!inp.is_open()){
            cerr << "Error reading file "<< filename <<endl;
            exit(1);
        }
        getstats(inp);
        inp.close ();
        ifstream inp1 (filename.c_str());
        read_stream (inp1, r, smap);
        inp1.close();
    }
}

void intervalmap::getstats (istream &inp ) {
	string line;
    io::println ("In getstats",2);
	while ( std::getline (inp, line)){
		char c = line[0];
		if (c=='#')
			continue;
        vector<string> toks;
        functions::tokenize(line, toks, " \t");

		string chr;
		double start;
		double end;
        
        chr = toks[0];
        start = atof(toks[1].c_str());
        end = atof(toks[2].c_str());

        if(givenchrom && chrom.compare (chr)!=0) 
            continue;

		if (chrmap.find(chr) == chrmap.end ()) {
            chrmap.insert (chr);
            intsperchr[chr] =  0;
			chrs.push_back(chr);
		}


        intsperchr[chr] = intsperchr[chr]+1;
    }
    io::println ("Out getstats",2);
    io::println ("chrsize = " + tostring (chrs.size()), 2);

	for (int i = 0 ; i < chrs.size(); i++) {
        string chr = chrs[i];
        cout << chr << "\t" << intsperchr[chr] << endl;
        intervals[chr] = vector<interval>(intsperchr[chr]);
    }
}

void intervalmap::read_stream (istream &inp, gsl_rng *r, snpmap *smap ) {
	string line;
	int id = 0;
    unordered_map <string,int> perchrindex;
    io::println ("read_stream",2);
    for (int i =  0 ;  i < chrs.size(); i++) 
        perchrindex[chrs[i]] =  0 ;

	while ( std::getline (inp, line)){
		char c = line[0];
		if (c=='#')
			continue;
//		istringstream ss (line);
        vector<string> toks;
        functions::tokenize(line, toks, " \t");

		string chr;
		double start;
		double end;
        
        chr = toks[0];
        start = atof(toks[1].c_str());
        end = atof(toks[2].c_str());


        if(givenchrom && chrom.compare (chr)!=0) 
            continue;

        string info =  "";
        for (int j = 3; j < toks.size();  j++)
            info += toks[j] + "\t";

        int n = perchrindex[chr];
		intervals[chr][n].set(id,chr,start,end,info);
        perchrindex[chr] = n+1;
		id ++;
	}

	for (int i = 0 ; i < chrs.size(); i++) {
		vector<snp>& tmp=snps[chrs[i]];
		sort (tmp.begin(), tmp.end());
		vector<interval>& tmpint = intervals[chrs[i]];
		sort(tmpint.begin(), tmpint.end());
		io::println ("Intervals on chr " + chrs[i] + "\t" + tostring(tmpint.size()), 2);
//		allsnps.insert (tmp.begin(), tmp.end());
	}
    /*
    vector<interval>& tmpint = intervals[chrs[0]];
    for (int i = 0 ;  i  < tmpint.size();  i++) { 
        cout.precision (0); cout.setf(ios::fixed); 
        cout << noshowpoint << tmpint[i].chr << "\t" << tmpint[i].start << "\t" << tmpint[i].end << endl;
    }*/

	this->r = r;
	this->smap = smap;
}


string intervalmap::stats () { 
	string s = "";
	
	unordered_map <string, vector<interval> >::iterator it;

	for (it = intervals.begin(); it != intervals.end (); it++) { 
		vector<interval> & tmpints = it->second;
		string chr = it->first;
		vector<snp> &snps = (smap->snps)[chr];
		for (int j = 0 ; j < tmpints.size(); j++){ 
			interval &i = tmpints[j];
			vector<int> &snpind = snpsperint[i];
			vector<int> &counts = countsperint[i];
			s += i.to_string() + "\twith " + tostring(counts.size()) + " intervals\n" ;
			for (int k = 0 ; k < counts.size(); k++){ 
				snp &snp = snps[snpind[k]];
				int c = counts[k];
				s = s +  snp.to_string() + "," + snp.stats() +",c=" + tostring(c) + "\n";
			}
			s += "\n";
		}
	}
	return s;
}

string intervalmap::to_string() const { 
	string s = "";
	for (int i = 0 ; i < chrs.size(); i++) {
		unordered_map<string, vector<interval> >::const_iterator it  = intervals.find(chrs[i]);
		const vector<interval>& tmp = it->second;
		for (int j = 0; j < tmp.size();j++) {
			s = s +  tmp[j].to_string ()  + "\n";
//            cout << j << "\t" << s <<endl;
        }
  //      cout <<i <<"\t" << s << endl;
	}
    return (s);
}

void intervalmap::add (const intervalmap &im ) {
    unordered_map<string,vector<interval> > imint = im.intervals;
	for (int i = 0 ; i < chrs.size(); i++) {
        if (imint.find (chrs[i])==imint.end())
                continue;
		unordered_map<string, vector<interval> >::const_iterator it  = intervals.find(chrs[i]);
		const vector<interval>& tmp1 = it->second;
        const vector<interval>& tmp2 = imint[chrs[i]];
        vector<interval> tmp;
        for (int j = 0, k = 0 ; ;) { 
            if (j < tmp1.size () && k < tmp2.size()) {
            if (tmp1[j] < tmp2[k]) {
                tmp.push_back(tmp1[j]);
                j++;
            } else {
                tmp.push_back(tmp2[k]);
                k++;
            }
            } else if (j < tmp1.size()) { 
                tmp.push_back(tmp1[j]);
                j++;
            } else if (k  <tmp2.size())  {
                tmp.push_back(tmp2[k]);
                k++;
            } else 
                break;
        }
        intervals[chrs[i]] = tmp;
    }	
    for ( int i  = 0 ; i < im.chrs.size() ;  i++) { 
        string chr = im.chrs[i];
        if (intervals.find(chr)==intervals.end()){
            intervals[chr] = vector<interval> (imint[chr]);
        }
    }
}



void intervalmap::intersect (const intervalmap &im ) {
    unordered_map<string,vector<interval> > imint = im.intervals;
	for (int i = 0 ; i < chrs.size(); i++) {
        if (imint.find (chrs[i])==imint.end()) {
                intervals.erase(chrs[i]);
                continue;
        }
		unordered_map<string, vector<interval> >::const_iterator it  = intervals.find(chrs[i]);
		const vector<interval>& tmp1 = it->second;
        const vector<interval>& tmp2 = imint[chrs[i]];
        vector<interval> tmp;
        for (int j = 0, k = 0 ; j < tmp1.size() && k < tmp2.size() ;) { 
            if (tmp1[j] && tmp2[k]) {
                interval tmp3 = tmp1[j]*tmp2[k];
                tmp.push_back(tmp3);
                if (tmp3.end==tmp1[j].end)
                    j++;
                else
                    k++;
            } else if (tmp1 < tmp2 && j < tmp1.size()) { 
                j++;
            } else if (tmp2 < tmp1 && k < tmp2.size()){
                k++;
            } else {
                break;
            }
        }
        intervals[chrs[i]]=tmp;
 
    }
}

void intervalmap::merge ( ) { 
	for (int l = 0 ; l < chrs.size(); l++) {
		unordered_map<string, vector<interval> >::const_iterator it  = intervals.find(chrs[l]);
		const vector<interval>& tmp = it->second;
        vector<interval> tmp1;
        for (int i =0 ; i < tmp.size () ;){ 
            int j = i + 1;
            interval tmp2 = tmp[i];
            while (j < tmp.size()){
                if ( tmp[i] && tmp[j]) { 
                    tmp2 = tmp2+tmp[j];
                    j++;
                } else { 
                    break;
                }
            }
            tmp1.push_back (tmp2);
            i = j;
        }
        intervals[chrs[l]] = tmp1;
    }
}
