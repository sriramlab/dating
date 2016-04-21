#include "std.h"
#include "functions.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "genotype.h"

ofstream plfs;
ofstream tlfs;
ofstream clfs;

class stratifysites { 
    public:
        string predfile;
        string snpfile;
        string genofile;
        string indfile;
        string indivfilterfile;
        bool givenindfile;
        bool givenindivfilterfile;
        string fillsnpfile;
        bool givenfillsnpfile;

        int verbose;
        int bins;
        int partition;
        double thresh1;
        double thresh2;
        double tlength;
        bool strict;
        bool usegenpos;
        
        bool givenprefix;
        string prefix;

        bool givennfreq;
        string nfreqfile;

        bool givenyfreq;
        string yfreqfile;

        bool givensep;
        string sep;

        bool giventrue;
        string truefile;


        int nind;
        int nsnps;
        int nchr;

        data *d;
        genotype *geno;
        packedgtype *predanc, *trueanc;
        unordered_map <string, double> *nfreqmap, *yfreqmap;

        stratifysites (int argc, char *argv[]) ;
        void run ();
};

stratifysites::stratifysites (int argc, char *argv[]) {
	static const char *optString = "vh";
	static const struct option longOpts[] = {
		{ "verbose", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
        {"genpos", no_argument, NULL, 'g'},
        {"predfile", required_argument, NULL, 'p'},
        {"snpfile", required_argument, NULL, 's'},
        {"bin", required_argument, NULL, 'b'},
        {"genofile", required_argument, NULL, 'g'},
        {"indfile", required_argument, NULL, 'i'},
        {"indivfilterfile", required_argument, NULL, 'i'},
        {"thresh1", required_argument, NULL, 't'},
        {"thresh2", required_argument, NULL, 't'},
        {"tlength", required_argument, NULL, 't'},
        {"strict", no_argument, NULL, 's'},
        {"fill", required_argument, NULL, 'f'},
        {"prefix", required_argument, NULL, 'g'},
        {"nfreq", required_argument, NULL, 'n'},
        {"yfreq", required_argument,  NULL, 'y'},
        {"sep", required_argument, NULL, 's'},
        {"truefile", required_argument, NULL,'t'},
		{ NULL, no_argument, NULL, 0 },
	};
	d = new data (argc, argv, optString, longOpts);
    d->get_string ("predfile", predfile, true);
    d->get_string ("snpfile", snpfile, true);
    d->get_int ("bin", bins, 10);
    d->get_string ("genofile", genofile, true);
    d->get_double ("thresh1", thresh1, 0.90);
    d->get_double ("thresh2", thresh2, 0.10);
    d->get_double ("tlength", tlength, 0);
    d->get_boolean ("strict", strict, false, false);

    givenindivfilterfile = d->get_string ("indivfilterfile", indivfilterfile, false);
    givenindfile = d->get_string ("indfile", indfile, false);

    givenfillsnpfile =  d->get_string ("fill", fillsnpfile, false);
    d->get_boolean ("genpos", usegenpos, false, false);
    givenprefix = d->get_string ("prefix", prefix, false );
    givennfreq = d->get_string ("nfreq", nfreqfile, false );
    givenyfreq = d->get_string ("yfreq", yfreqfile, false );
    givensep = d->get_string ("sep", sep, false );
    giventrue = d->get_string ("truefile", truefile, false );

    if (givenindivfilterfile && givenindfile)
        geno = new genotype ( snpfile, indfile, genofile, "eigenstrat", false, indivfilterfile, "", true);
    else
        geno = new genotype ( snpfile, genofile, false, false, true);

    nind = geno->nind;
    nsnps = geno->nsnp;
    nchr = geno->nchr;

    predanc = new packedgtype (nind, nsnps, geno->cumnumsnps, false, false, NULL, NULL, "prob");
    predanc->read_eigenstrat (predfile, NULL, NULL);
    if (io::debug>=2) { 
        for ( int i =  0 ; i < geno->nchr; i++) { 
            for (int j = 441876 ; j < nsnps && j < 441976 ; j++) { 
                for (int l = 0 ; l < nind ; l++) { 
                    double g = (*predanc).get(i,j,l);
                    cout << g << "\t";
                }
                cout <<endl;
            }
        }
    }

    if( givennfreq) {
        unordered_map<string, string> *tmpmap1 = new unordered_map<string, string> ();
        nfreqmap = new unordered_map<string, double> ();
        fileio::read_map (nfreqfile,tmpmap1,0,1);
        for (unordered_map<string,string>::iterator i = tmpmap1->begin(); i != tmpmap1->end(); i++) { 
            (*nfreqmap)[i->first] = atof ((i->second).c_str());
        }
    }

    if( givenyfreq) {
        unordered_map<string, string> *tmpmap1 = new unordered_map<string, string> ();
        yfreqmap = new unordered_map<string, double> ();
        fileio::read_map (yfreqfile,tmpmap1,0,1);
        for (unordered_map<string,string>::iterator i = tmpmap1->begin(); i != tmpmap1->end(); i++) { 
            (*yfreqmap)[i->first] = atof ((i->second).c_str());
        }
    }

    if (giventrue) { 
        trueanc = new packedgtype (nind, nsnps, geno->cumnumsnps, false);
        trueanc->read_eigenstrat (truefile, NULL, NULL);
    }

}

void stratifysites::run ()  {
    cout << "# Parameters\n";
    cout << "# thresh1=" << thresh1 << "\n";
    cout << "# thresh2=" <<thresh2 << "\n";
    cout << "# tlength="<<tlength << "h\n";
    cout << "# genpos="<<usegenpos << "\n";
    cout << "# strict="<<strict << "\n";
    cout << "# bin="<<bins << "\n";
    cout << "# predfile="<<predfile << "\n";
    cout << "# snpfile="<<snpfile << "\n";
    cout << "# genofile="<<genofile << "\n";
    cout <<  "## SNP id\n";
    cout <<  "## Chromosome \n";
    cout <<  "## Genetic position\n";
    cout <<  "## Physical position\n";
    if (givennfreq)
        cout <<  "## Neandertal DAF\n" ;
    if (givenyfreq)
        cout <<  "## YRI DAF\n" ;
    cout <<  "## Number of derived alleles\n";
    cout <<  "## All derived alleles on neandertal haplotype ? \n";
    cout <<  "## Some derived alleles on neandertal and some on human ?\n";
    cout <<  "## Average posterior probability of N on derived alleles \n";
    cout <<  "## Average posterior probability of N \n";
    cout <<  "## Derived allele frequency\n";
    cout <<  "## Neandertal haplotype frequency\n";
    cout <<  "## DAF on N haplotypes\n";
    cout <<  "## Number of N haplotypes\n";
    cout <<  "## DAF on MH haplotypes\n";
    cout <<  "## Number of MH haplotypes\n";



    if (givenprefix) {
        plfs.open (prefix+".haplotypes");
        tlfs.open (prefix+".contigs");
        clfs.open ( prefix +".lanc");

        tlfs.setf(ios::fixed,ios::floatfield); 
        tlfs.precision(4);
        tlfs <<  "## Chromosome\n";
        tlfs << "## Contig start\n";
        tlfs << "## Contig end\n";
        tlfs << "## Average neandertal ancestry\n";
        tlfs << "## Number of SNPs\n";
        
        plfs.setf(ios::fixed,ios::floatfield); 
        plfs.precision(4);
        plfs << "## Chromosome " <<endl;
        plfs << "## Individual " <<endl;
        plfs << "## Start pos " <<endl;
        plfs << "## End pos " <<endl;
        plfs << "## Physical length " <<endl;
        plfs << "## Genetic length " <<endl;
        plfs << "## Average predicted probability " <<endl;
        plfs << "## Average true probability " <<endl;
        plfs << "## Number of SNPs " <<endl;

        cout.setf(ios::fixed,ios::floatfield); 
        cout.precision(4);
        clfs.setf(ios::fixed,ios::floatfield); 
        clfs.precision(4);
        clfs << "## SNP id " <<endl;
        clfs << "##  Maximum physical length of overlapping N haplotype " << endl;
        clfs << "##  Maximum genetic length of overlapping N haplotype " << endl;
        clfs << "## Average physical length of overlapping N haplotype " << endl;
        clfs << "## Average genetic length of overlapping N haplotype " << endl;
    }

    for (int c = 0 ; c < nchr; c++) { 
        vector<snp> snps = geno->snps[geno->chrs[c]];
        vector<int> tmppred (nind, 0);
        vector<double> reg (nind,0);
        vector<double> reg2 (nind,0);
        vector<double> st1 (nind,0);
        vector<double> st2 (nind,0);
        vector<int> startpos (nind, 0);
        
        vector<double> maxlengths(snps.size(), 0 );
        vector<double> avglengths(snps.size(), 0 );
        vector<double> maxglengths(snps.size(), 0 );
        vector<double> avgglengths(snps.size(), 0 );
        vector<int> denomlengths(snps.size(), 0 );

        unordered_map <string, double> plmap ;
        unordered_map <string, double> glmap ;

        int incontig = 0 ;
        int contigsnps  = 0;
        string contigchr;
        double contigstart = 0;
        double  contigend = 0 ;
        double contigl = 0 ;
        double contigavganc =  0;

        double avgpredanc  = 0;

        for (int i = 0 ; i < snps.size() ; i++) { 
            for (int j = 0 ; j < nind; j++) { 
                double p = (*predanc).get(c,i,j);
                tmppred [j] = (p>thresh1)?1:0;
                avgpredanc += p;
            }
            string id = snps[i].id;
            string chr = snps[i].chr;

            if (nind>0)
                avgpredanc /= nind ;

            if (i==0) { 
                for (int j = 0  ; j < nind ; j++) { 
                    reg[j]  =reg2[j] = 0;
                    st1[j] = st2[j] = tmppred[j];
                    startpos[j] = 0;
                }
            } else { 
                int count1 = 0 ; 
                int count2 = 0;
                for (int j = 0 ; j < nind; j++) { 
                    int flag1 = 0;
                    int flag2 = 0;
                    int len1 = 0 ;
                    double len2 = 0;
                    double avgconf = 0 ;
                    double trueconf = 0 ;
                    int denom = 0 ;
                    int start = 0 ;
                    int end = 0 ;

                    if (st1[j] == 0 && tmppred[j] > 0) {
                        st1[j] = 1;
                        startpos[j] = i;
                    } else if (st1[j]==1 && tmppred[j]==0) {
                        len1 = reg[j];
                        for (int k = i - 1; k >= 0 && k >= startpos[j] ; k--) {
                            if ( maxlengths[k] < len1) {
                                maxlengths[k] = len1;
                            }
                            avglengths[k] += len1;
                            denomlengths[k]++;

                            start = snps[k].physpos;
                            avgconf +=  (*predanc).get (c,k,j);
                            if (giventrue)
                                trueconf +=  (*trueanc)(c,k,j);
                            denom ++;

                            string key =  chr + ":" + tostring(k) +":"+ tostring(j);
                            plmap[key] = len1;

                        }
                        avgconf = (denom>0)?avgconf/denom:0;
                        trueconf = (denom>0)?trueconf/denom:0;
                        reg[j] = 0;
                        st1[j] = 0;
                        flag1 = 1;
                        end = snps[i-1].physpos;
                        if (len1 > 0) count1++;
                    } else if (st1[j]==1 && tmppred[j] > 0 ) {
                        reg[j] +=  tmppred[j] *(snps[i].physpos - snps[i-1].physpos);
                        count2 ++;
                    } else {
                    }


                    if (st2[j] == 0 && tmppred[j] > 0)
                        st2[j] = 1;
                    else if (st2[j]==1 && tmppred[j]==0) {
                        len2 = reg2[j];
                        for (int k = i - 1; k >= 0 && k >= startpos[j] ; k--) {
                            if (maxglengths[k] < len2) {
                                maxglengths[k] = len2;
                            }
                            avgglengths[k] += len2;
                            string key =  chr + ":" + tostring(k) +":"+ tostring(j);
                            glmap[key] = len2;
                        }
                        reg2[j] = 0; 
                        st2[j] = 0 ;
                        flag2 = 1;

                    } else if (st2[j]==1 && tmppred[j] > 0 ) {
                        reg2[j] +=  tmppred[j] *(snps[i].genpos - snps[i-1].genpos);
                    } else {
                    }

                    if (flag1 && flag2) {
                        if (givenprefix) {

                            len2 *= 100;
                            plfs << geno->chrs[c] << "\t" << j << "\t" << start << "\t" << end <<"\t" << len1 << "\t" << len2;
                            plfs <<"\t"<< avgconf << "\t" << trueconf << "\t" << denom << endl;
                        }
                    }
                }
                if (count2 > 0  ) { 
                    if (incontig) { 
                        contigsnps ++;
                        contigavganc += avgpredanc;
                    } else { 
                        contigchr = geno->chrs[c];
                        contigstart = snps[i-1].physpos;
                        incontig = 1;
                    }
                    contigl += (snps[i].physpos - snps[i-1].physpos);
                } else if (count2 == 0 && count1 >  0) { 
                    contigend = snps[i].physpos;
                    if (contigsnps > 0) 
                        contigavganc /= contigsnps;
                    long cs = contigstart;
                    long ce = contigend;
                    if (givenprefix)
                        tlfs << contigchr << "\t" << cs << "\t" << ce << "\t" << contigavganc << "\t" << contigsnps << endl;
                    contigl = 0 ;
                    incontig = 0 ;
                    contigsnps = 0 ;
                    contigavganc = 0 ;
                }
            }
        }

        for (int i = 0 ; i < snps.size() ; i++) { 
            string id = snps[i].id;
            string chr = snps[i].chr;
            double pos = snps[i].physpos;
            double gpos = snps[i].genpos;

            int n1 = 0;
            int n2 = 0;
            double n3 = 0;
            double n4 = 0;
            double g = 0 ;
            double f = 0 ;
            double f2 = 0 ; 
            int denom2 = 0 ;
            double f1 = 0 ;
            int denom1 = 0 ;
            int d = 0 ;


            double maxl = 0 ;
            double maxgl  = 0;
            for (int j = 0 ; j < nind ; j++) { 
                string key = chr + ":" + tostring(i) + ":" + tostring(j);
                double panc = (*predanc).get(c,i,j);
                if (io::debug >= 1)
                    cout << "panc ( " << c <<","<<i<<","<<j<<") = " << panc << endl;
                int predn = 0 ;
                int predh = 0 ;
                double predl = 0 ;
                if (usegenpos) {
                    if ( glmap.find(key) != glmap.end())
                        predl = glmap[key] * 100;
                    if (maxgl < predl)
                        maxgl = predl;
                } else { 
                    if ( plmap.find(key) != plmap.end())
                        predl = plmap[key];
                    if (maxl < predl) 
                        maxl  = predl;
                }

                if ( panc >= thresh1 && predl >= tlength)
                    predn = 1;
                if (panc <= thresh2)
                    predh = 1;

                int allele = (*geno)(c,i,j);
                f += allele;
                if (allele==1) {
                    d ++;
                    n1 += predn;
                    n2 += predh;
                    n3 += panc;
                }
                n4 += panc;

                if (predn) {
                    f1 += allele;
                    denom1++;
                }
                if (predh) { 
                    f2 += allele;
                    denom2 ++;
                }

            }
            f /= nind;
            g /= nind;
            f1 = denom1>0?(f1/denom1):0;
            f2 = denom2>0?(f2/denom2):0;

            if (d>0) { 
                n1/=d; n2/=d; n3/=d; 
            }
            if (nind > 0)
                n4 /= nind;
            int m1 = (n1>=1)?1:0;
            int m2 = ( (n1>0)&&(n2>0))?1:0;

            cout << id << "\t" << chr << "\t" << gpos*100 << "\t" << snps[i].getphyspos ();
            if (givennfreq )
                cout << "\t" << (*nfreqmap)[id];
            if (givenyfreq )
                cout << "\t" << (*yfreqmap)[id];
            cout << "\t" << d << "\t" << m1 << "\t"<< m2 << "\t" << n3 << "\t" << n4;
            cout << "\t" << f << "\t" << g <<"\t" << f1 <<"\t" << denom1 <<"\t" << f2 <<"\t" << denom2;
            cout << endl;

            if (denomlengths[i]>0) { 
                avglengths[i] /= denomlengths[i];
                avgglengths[i] /= denomlengths[i];
            }

            if (givenprefix) {
                maxglengths[i] *= 100;
                avgglengths[i] *= 100;
                long ml = maxlengths[i];
                long al = avglengths[i];
                clfs << id << "\t" << ml  << "\t" << maxglengths[i] << "\t" << al << "\t" << avgglengths[i] << "\t" <<denomlengths[i] <<endl;
            }

        }


    }
    if (givenprefix) {
        plfs.close ();
        tlfs.close ();
        clfs.close ();
    }

}

int main (int argc, char* argv[]) {
	stratifysites s = stratifysites(argc,argv);
    s.run ();
}
