#include "intervalmap.h"
#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "vectorfn.h"
#include "snpmap.h"
#include "convertf.h"



class intervalstats  {

    public:
        intervalmap *imap;
        intervalstats (int argc, char *argv[]); 
        
        // Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;


};

intervalstats::intervalstats (int argc, char *argv[])  {
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


    p->get_int ("debug", io::debug, 0);
    string intervalfile;
    p->get_string ("intervals",intervalfile, true);
    imap = new intervalmap (intervalfile, NULL, NULL);



    string freqfile;
    p->get_string ("freqfile", freqfile, true);

    string outfile;
    bool givenoutfile = false;
    givenoutfile = p->get_string ("outfile", outfile, false);

    unordered_map<string, vector<interval> > &intervals = imap->intervals;
    ifstream ifs (freqfile) ;

    string line;
    string prevchr = "";
    string chr ; 
    double pos ;
    double stat;
    double prevpos ; 
    double prevstat;
    int linenum =  0;
    int intervalind = 0;
    vector<interval> *intervalsonchr;
    vector<int> activeset;
    vector<double> startpos;
    vector<double> stats;
    vector<int> numsnps;
    ofstream outfs;
    if (givenoutfile)
        outfs.open (outfile.c_str());
    while (std::getline (ifs, line)) {
        vector<string> toks;
        functions::tokenize (line, toks, " \t");
        chr = toks[0];
        pos = atof(toks[1].c_str());
        stat = atof(toks[2].c_str());
        if (io::debug >= 1)  { 
            cout << chr << "\t" << pos << endl;
        }
        if (linenum==0 || chr.compare(prevchr)!=0) {

            if (linenum >0 ) { 
                for (int i = 0  ; i < activeset.size(); i++) {
                    interval &tmpint = (*intervalsonchr)[activeset[i]];
                    double x = stats[i]/(prevpos - startpos[i]);
                    if (givenoutfile) 
                        outfs << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << x <<"\t" << numsnps[i]  << endl;
                    else
                        cout << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << x << "\t" << numsnps[i]  << endl;
                }
                for (int i = intervalind; i < intervalsonchr->size(); i++)  {
                    interval &tmpint = (*intervalsonchr)[i];
                    if (givenoutfile) 
                        outfs << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t0\t0"  << endl;
                    else
                        cout << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t0\t0"  << endl;
                }
                activeset.clear ();
                startpos.clear();
                stats.clear ();
                numsnps.clear ();
            }

            intervalsonchr = &(intervals[chr]);
            for (int i =  0 ;i < intervalsonchr->size() ; i++) { 
                interval &tmpint = (*intervalsonchr)[i];
                if ( tmpint.start > pos) {
                    intervalind =  i ;
                    break;
                }
                if (tmpint.end <= pos) { 
                    if (givenoutfile) 
                        outfs << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t0\t0"  << endl;
                    else
                        cout << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t0\t0"  << endl;
                } else {
                    activeset.push_back (i);
                    startpos.push_back (pos);
                    stats.push_back (0);
                    numsnps.push_back (1);
                }
            }

            prevpos = pos ;
            prevstat = stat;
            prevchr = chr;
            linenum++;
            continue;

        } else {
        }


        if (io::debug >= 1)  {
            if (activeset.size() > 0) {
                cout << "activeset1 = " ;
                vectorfn::printvector (activeset);
            }
        }

        vector<int> tmpactiveset;
        vector<double> tmpstats;
        vector<double> tmpstart;
        vector<int> tmpnumsnps;
        for (int i =  0 ; i < activeset.size (); i++) {
            interval &tmpint = (*intervalsonchr)[activeset[i]];
            if (tmpint.end <= pos ) { 
                stats[i] += 0.5*(tmpint.end-prevpos)*(prevstat+stat);
                stats[i] /= (tmpint.end - startpos[i]);
                if (givenoutfile)
                    outfs << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << stats[i]  <<"\t"<<numsnps[i]<< endl;
                else
                    cout << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << stats[i]  <<"\t" << numsnps[i]<< endl;

            } else {
                stats[i] += 0.5*(pos-prevpos)*(prevstat+stat);
                tmpactiveset.push_back (activeset[i]);
                tmpstats.push_back (stats[i]);
                tmpstart.push_back (startpos[i]);
                tmpnumsnps.push_back (numsnps[i]+1);
            }
        }

        activeset.clear (); 
        stats.clear ();
        startpos.clear ();
        numsnps.clear (); 
        if (tmpactiveset.size() >  0)  {
            activeset = tmpactiveset; 
            stats = tmpstats;
            startpos  = tmpstart;
            numsnps  = tmpnumsnps;
        }
        if (io::debug >= 1)  {
            if (activeset.size() > 0) {
                cout << "activeset2 = " ;
                vectorfn::printvector (activeset);
                cout << "intervals = " << intervalsonchr->size() << endl;
            }
        }


        for (int i =  intervalind ;i < intervalsonchr->size() ; i++) { 
            interval &tmpint = (*intervalsonchr)[i];
            if ( tmpint.start > pos) {
                intervalind =  i;
                break;
            }
            double statatstart = (tmpint.start-prevpos)*prevstat  + (pos-tmpint.start)*stat;
            statatstart /= (pos-prevpos);

            if (tmpint.end <= pos) {
                double statatend = (tmpint.end-prevpos)*prevstat  + (pos-tmpint.end)*stat;
                statatend /= (pos-prevpos);
                double tmp = 0.5 * ( tmpint.end-tmpint.start)*(statatstart + statatend);
                tmp /= (tmpint.end  - tmpint.start);
                if (givenoutfile)
                    outfs << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << tmp  << "\t0" << endl;
                else
                    cout << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << tmp  << "\t0"<<endl;

            } else {
                double tmp = 0.5*(pos-tmpint.start)*(statatstart+stat);
                activeset.push_back (i);
                startpos.push_back (tmpint.start);
                stats.push_back (tmp);
                numsnps.push_back (1);
            }
            intervalind ++;
        }

        prevchr = chr;
        prevstat = stat;
        prevpos =  pos;
        linenum++;
        if (io::debug >= 1)
            cout << "intervalind = " << intervalind << "\t intervalsonchr = " << intervalsonchr->size() << endl;
    }
    for (int i = 0  ; i < activeset.size(); i++) {
        interval &tmpint = (*intervalsonchr)[activeset[i]];
        double x = stats[i]/(prevpos - startpos[i]);
        if (givenoutfile)
            outfs << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << x <<"\t" << numsnps[i] << endl;
        else
            cout << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t" << x  <<"\t" << numsnps[i] << endl;
    }
    for (int i = intervalind; i < intervalsonchr->size(); i++)  {
        interval &tmpint = (*intervalsonchr)[i];
        if (givenoutfile)
            outfs << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t0\t0"  << endl;
        else
            cout << noshowpoint << tmpint.chr << "\t" << snp::tophyspos (tmpint.start) << "\t" << snp::tophyspos(tmpint.end) << "\t" << tmpint.info << "\t0\t0"  << endl;
    }

}



int main (int argc , char *argv[]) { 
    intervalstats f = intervalstats (argc, argv);
}
