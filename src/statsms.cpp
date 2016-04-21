#include "data.h"
#include "io.h"
#include "functions.h"
#include "stringfn.h"
#include "mathfn.h"
#include "statsms.h"
#include "tajd.h"
#include "afs.h"

ofstream outfs;
ofstream disfs;
ofstream divfs;

statsms::statsms (string configfile) { 
	string s;
	d = new data (configfile);
	d->get_int ("debug", io::debug, 0);
	d->get_string ("param", paramfile, true);
	d->get_string ("ms",msfile,true);
	param = new data (paramfile);
	param->get_int ("length", length, 1000000, false);
	param->get_double ("rho", rho, 1.3e-8, false);
	param->get_int ("reps", reps, 100, true);
	param->get_string ("samples", s, true);
    string s1 = "ooa,afr,n,d,c";
    param->get_string ("samplenames", s1, false);

	io::println ("samples=  " + s ,2);
	vector<string> toks; 
	functions::tokenize (s,  toks, ",");
	for (int i = 0 ; i < toks.size(); i++)
		samples.push_back (atoi (toks[i].c_str()));
	
    toks.resize(0);
    functions::tokenize(s1,toks,",");
    samplenames.resize (toks.size());
    for (int i = 0 ; i <toks.size(); i++) { 
        samplemap[toks[i]]=i;
        samplenames[i] = toks[i];
    }
    ooa = samplemap["ooa"];
    afr = samplemap["afr"];
    den = samplemap["d"];
    nea = samplemap["n"];
    chimp = samplemap["c"];


	nsamples = samples.size();
	for (int i  = 0 ; i < samples.size();i++)
		io::print( tostring(samples[i])+ ",",2);
	io::print("\n",2);


	d->get_double ("maxdis", maxdis, 1e-2);
	d->get_double ("binsize", binsize,1e-5);
	double x = maxdis/binsize;
	nbins = mathfn::round (x);
	io::println ("nbins="+tostring (nbins),2);

	d->get_int ("derived",derived, 1,false);
	d->get_int ("poly",poly, 1,false);
	d->get_double ("freq", fcutoff, 0.1);
    d->get_int ("replicates", replicates, 100);
    if (replicates > reps)
        replicates = reps;
    chrsperrep = reps/replicates;

//	d->get_int ("ascertain", ascertain, 2, false);
	d->get_int ("uniform", uniform, 1, false);
	d->get_boolean ("fast", fast, true, false);
	d->get_boolean ("asn", asn, false, false);
	d->get_boolean ("kg", kg, false, false);
	d->get_int ("seed",seed,1,false);
	d->get_int ("usegenetic",usegenetic,0,false);
	d->get_boolean ("mshot", mshot,false,false);
	d->get_double ("alpha", alpha, 1e3, false);
	d->get_int ("start", start, 0, false);
	d->get_int ("corr", corr, 0, false);
    d->get_boolean ("hetascertain", hetascertain, true, false);
    ascertain = "outgroup";
    d->get_string ("ascertain", ascertain,  false);
	printsnplist =  d->get_string ("snpfile", snpfile, false);

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, seed);


	out = d->get_string ("outfile", outfile, false);
	if (out) {
	       outfs.open (outfile.c_str());	
	       string disfile = outfile + ".dis";
	       disfs.open (disfile.c_str());
	}
    if (io::debug >=1)
    	d->print_parameters();

    divfs.open ("div.txt");
}


void statsms::readms (string filename){  
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
	int chrindex = 0;
	int snpindex = 0;
    int nsnp = 0;
    double totalw =  0;
	ofstream snpfs (snpfile.c_str());

	
    int totalsnps =  0;
    vector<int> nsegsites (nsamples,0);
    int totalnind = 0 ;
    vector<int> nind  (nsamples,0);
    for (int i = 0  ; i < nsamples; i++) {
        nind[i] = samples[i];
        totalnind += nind[i];
    }

    afs jafs ( nind, samplenames);

    char ** snpmat = new char*[totalnind];
    double **freq = new double*[nsamples];
    double **td = new double *[nsamples];
    double **thetaw = new double *[nsamples];
    for (int i  = 0 ; i < nsamples; i++) { 
        td[i] = new double[reps];
        thetaw[i]  = new double[reps];
    }
    double abba =  0 ;
    double baba = 0 ;
    double abba2 =  0 ;
    double baba2 = 0 ;
    double fstnum  = 0;
    double fstdenom =  0;
    double l2yn = 0 ;
    double l2en =  0;
    double l2denom  =  0;

    vector<double> pop1nd10 (nind[0]+1,0);
    vector<double> pop1nd01 (nind[0]+1,0);
    vector<double> pop2nd10 (nind[1]+1,0);
    vector<double> pop2nd01 (nind[1]+1,0);
    vector<double> pop1nd00 (nind[0]+1,0);
    vector<double> pop2nd00 (nind[1]+1,0);
    vector<double> pop1nd11 (nind[0]+1,0);
    vector<double> pop2nd11 (nind[1]+1,0);
    double divan = 0 ;
    double divad = 0 ;

    vector< vector<double> > jpop1nd00 (replicates, vector<double> (nind[0]+1,0) );
    vector< vector<double> > jpop1nd01 (replicates, vector<double> (nind[0]+1,0));
    vector< vector<double> > jpop1nd10 (replicates, vector<double> (nind[0]+1,0));
    vector< vector<double> > jpop1nd11 (replicates, vector<double> (nind[0]+1,0));
    vector< vector<double> > jpop2nd00 (replicates, vector<double> (nind[1]+1,0));
    vector< vector<double> > jpop2nd01 (replicates, vector<double> (nind[1]+1,0));
    vector< vector<double> > jpop2nd10 (replicates, vector<double> (nind[1]+1,0));
    vector< vector<double> > jpop2nd11 (replicates, vector<double> (nind[1]+1,0));
    vector<double> jdenom (replicates,0);
    vector<double> jdivan (replicates,0);
    vector<double> jdivad (replicates,0);
    vector<double> jdivr (replicates, 0 );




    vector<double> chrabba = vector<double>(reps,0);
    vector<double> chrbaba = vector<double>(reps,0);
    vector<double> chrsnps = vector<double>(reps,0);
    vector<double> chrfstnum = vector<double> (reps,0);
    vector<double> chrfstdenom = vector<double> (reps,0);
    vector<double> chrw = vector<double> (reps,0);
    int totalnsnps = 0;

	int rep = 0 ;
	int index = 0;
    int repindex =  -1;
	string id;
    double ndiv1 = 0 ; double ddiv1 = 0 ;

	while ( std::getline (inp, line)){
		linenum ++;
		io::println ("line = " + line, 2);
		char c = line[0];
		if (c=='#')
			continue;
		
		if (line.empty())
			continue;

        if (line.find("(")!=std::string::npos && line.find(")")!=std::string::npos) {
            ndiv1 = 0 ;ddiv1 = 0 ;
            cout << rep << "\tline = " << line << endl;
            vector<string> toks (0);
            functions::tokenize (line.c_str(), toks, "(),:");
            for (int i  = 0 ; i  < toks.size() ;  i++){
                cout << toks[i] << ";";
            }
            cout << endl;

            for (int i  = 0 ; i  < toks.size() ;  i++){
                if (toks[i].find (".")==std::string::npos){
                    int ind = atoi (toks[i].c_str());
                    if (ind==totalnind-2) {
                        ndiv1 = atof (toks[i+1].c_str());
                    }
                    if (ind==totalnind-1) {
                        ddiv1 = atof (toks[i+1].c_str());
                    }
                }
            }
        }

		if ( line.find ("segsites:")==0) { 
			string l1 = line.substr(9,line.length()-9);
			nsnp = atoi (l1.c_str());
            for (int i = 0 ; i < totalnind ; i++) {
                snpmat[i] = new char[nsnp];
            }

            for (int i = 0 ; i < nsamples; i++){
                freq[i] = new double[nsnp];
            }

			io::println ("SNPs = " + tostring(nsnp),2);
		}

		if ( line.find ("positions:")==0) { 
			vector<double> genmap ;
			vector<double> physmap;
			vector<string> snpids;
			vector<vector<int> > geno;
			vector<vector<int> > haplo;
			chrindex ++;
			
			string l1 = line.substr(10,line.length()-10);
			//cout << "l1 = " << l1 << endl;
			//
			//cout << " rep = " << rep  << endl;
            if (rep%chrsperrep == 0 && repindex < replicates ) 
                repindex ++;
			io::println  ("rep = " + tostring(rep),2);
			io::println  ("repindex = " + tostring(repindex),2);

			rep++;
			istringstream ss (line.substr(10,line.length()-10));
			int count = 0 ;
			while (!ss.eof()){
				double pos;
				ss>>pos;
			//	cout << pos << endl;
				double ppos = pos * length;
				double gpos = ppos * rho;
				if (count >= nsnp ) 
					break;
				count++;

				genmap.push_back (gpos);
				physmap.push_back (ppos);
				id = "snp" + tostring(index);
				snpids.push_back (id);
				index++;
			}
            for (int i = 0 ; i  < totalnind ; i ++){
				string s1, s2;
				std::getline (inp, s1);
				for (int j = 0 ; j < nsnp ; j++) {
                    snpmat[i][j] = s1[j];
                }
            }
            vector<int> nsegsites(nsamples,0);
            vector<double> pi(nsamples, 0 );
            vector<double> h(nsamples,0);
            vector<double> th(nsamples,0);
            for (int i = 0, t =  0 ; i < nsamples; i++){ 
                for (int j =  0 ; j < nsnp; j++){
                    freq[i][j] = frequency('1', j, samples[i], snpmat+t)/((double)samples[i]);
                    nsegsites[i] += (freq[i][j] > 0 && freq[i][j] < 1.0);
                }
                t += samples[i];
            }
            
            int nout = 0 ; int dout = 0 ;
            int nin = 0; int din = 0 ;
            int ndiv =  0 ; int ddiv = 0 ;
            for (int j = 0 ; j < nsnp ; j++)  {
                double a = (samples[afr]*freq[afr][j] + samples[den]*freq[den][j])/(samples[afr]+samples[den]);
                double b = (samples[afr]*freq[afr][j] + samples[nea]*freq[nea][j])/(samples[afr]+samples[nea]);
                if (freq[chimp][j]==0 && freq[nea][j]==0 && freq[afr][j]==afr && freq[den][j]==afr){
                    nout ++;
                }
                if (freq[chimp][j]==0 && freq[den][j]==0 && freq[afr][j]==afr && freq[nea][j]==afr){
                    dout ++;
                }
                if (freq[chimp][j]==0 && freq[nea][j]==afr && a > 0  && a < 1){
                    nin ++;
                }
                if (freq[chimp][j]==0 && freq[den][j]==afr && b > 0 && b < 1 ){
                    din ++;
                }
                if (freq[chimp][j]==0 && freq[afr][j]==0 && freq[den][j]==0 && freq[nea][j]==afr){
                    ndiv ++;
                }
                if (freq[chimp][j]==0 && freq[afr][j]==0 && freq[nea][j]==0 && freq[den][j]==afr){
                    ddiv ++;
                }
            }
            if (nout > 0 && nin > 0) {
                cerr << "Inconsistent " << endl;
            }
            if (dout > 0 && din > 0) {
                cerr << "Inconsistent " << endl;
            }
            bool noutflag = false; bool doutflag = false;
            noutflag = (nout>0 && nin==0);
            doutflag = (dout>0 && din==0);

            divfs << noutflag << "\t" << doutflag << "\t" << "\t" << ndiv << "\t" << ddiv << "\t" << ndiv1 << "\t" << ddiv1 << endl;

            chrsnps[rep-1] = nsnp;
            totalnsnps += nsnp;
            vector<int> config (nsamples);
            for (int j = 0 ; j < nsnp ; j++)  {

                if (nsamples >= 5 ) {
                    if (ascertain.compare ("african")==0) { 
                        if (freq[2][j] == 0  && freq[3][j] == 0 && freq[4][j] ==0 && freq[1][j]==0)
                            continue;
                    }
                }

                if (freq[2][j]>=1){
                    abba += freq[1][j] * ( 1-freq[0][j]) * freq[2][j];
                    baba += freq[0][j] * ( 1-freq[1][j]) * freq[2][j];
                    chrabba[rep-1] += freq[1][j] * ( 1-freq[0][j]) * freq[2][j];
                    chrbaba [rep-1] += freq[0][j] * ( 1-freq[1][j]) * freq[2][j];
                }
                

                for (int i = 0  ;  i < nsamples; i++) {
                    int ind = round(freq[i][j]*samples[i]);
                    config[i] = ind;
                }
                jafs.increment (config);
                    
                if (nsamples >= 4) { 
                    abba2 += freq[den][j]*freq[afr][j]*(1-freq[ooa][j]);
                    baba2 += freq[den][j]*(1-freq[afr][j])*freq[ooa][j];
                    int ind = round(freq[ooa][j]*samples[ooa]);
                    pop1nd00[ind] += (1-freq[nea][j])*(1-freq[den][j]);
                    pop1nd01[ind] += (1-freq[nea][j])*freq[den][j];
                    pop1nd10[ind] += freq[nea][j]*(1-freq[den][j]);
                    pop1nd11[ind] += freq[nea][j]*freq[den][j];
                    
                    jpop1nd00[repindex][ind] += (1-freq[nea][j])*(1-freq[den][j]);
                    jpop1nd01[repindex][ind] += (1-freq[nea][j])*freq[den][j];
                    jpop1nd10[repindex][ind] += freq[nea][j]*(1-freq[den][j]);
                    jpop1nd11[repindex][ind] += freq[nea][j]*freq[den][j];



                    ind = round(freq[afr][j]*samples[afr]);
                    pop2nd00[ind] += (1-freq[nea][j])*(1-freq[den][j]);
                    pop2nd01[ind] += (1-freq[nea][j])*freq[den][j];
                    pop2nd10[ind] += freq[nea][j]*(1-freq[den][j]);
                    pop2nd11[ind] += freq[nea][j]*freq[den][j];

                    jpop2nd00[repindex][ind] += (1-freq[nea][j])*(1-freq[den][j]);
                    jpop2nd01[repindex][ind] += (1-freq[nea][j])*freq[den][j];
                    jpop2nd10[repindex][ind] += freq[nea][j]*(1-freq[den][j]);
                    jpop2nd11[repindex][ind] += freq[nea][j]*freq[den][j];

                    jdenom[repindex]=jdenom[repindex]+1;

                    divan +=  (freq[afr][j]*(1 - freq[chimp][j]) *( 1-freq[nea][j]));
                    divad +=  (freq[afr][j]*(1 - freq[chimp][j]) *( 1-freq[den][j]));
                    jdivan[repindex] +=  (freq[afr][j]*(1 - freq[chimp][j]) *( 1-freq[nea][j]));
                    jdivad[repindex] +=  (freq[afr][j]*(1 - freq[chimp][j]) *( 1-freq[den][j]));
                }


                if (freq[0][j]<=0 &&& freq[1][j]>0 && freq[1][j]<1)  {
                    l2en +=pow (freq[0][j] - freq[2][j],2);
                    l2yn +=pow (freq[1][j] - freq[2][j],2);
                    l2denom ++;
                }

                int n0 = samples[0];
                int n1 = samples[1];
                int a0 = freq[0][j] * n0;
                int a1 = freq[1][j] * n1 ;
                double w = 1;
                double f0 = freq[0][j];
                double f1 = freq[1][j];
                if (ascertain.compare("hetinb")==0) {  
                    a1 --; // Account for ascertainment
                    n1 -= 2;
                    f1 = ((double)a1)/n1;
                    w = 2*freq[1][j]*(1-freq[1][j]);
                } else if (ascertain.compare("outgroup")==0) { 
                    if (freq[2][j]>0)
                       w = 1 ;
                    else 
                       w = 0; 
                } else if (ascertain.compare ("all")==0) {
                    w = 1;
                } else if (ascertain.compare ("hetinout") == 0) { 
                    w = freq[2][j]*(1-freq[2][j]);
                } else if (ascertain.compare ("polyinb")==0) {
                    if (freq[1][j]>0 && freq[1][j] < 1 ) w = 1;
                    else w = 0;
                } else if (ascertain.compare ("polyinanb")==0){ 
                    if (freq[1][j]>0 && freq[1][j] < 1 && freq[0][j] > 0 && freq[0][j] < 1) w = 1;
                    else w = 0;
                } else if (ascertain.compare ("polyinaorb")==0) { 
                    if ((freq[1][j]>0 && freq[1][j] < 1)||(freq[0][j] > 0 && freq[0][j] < 1)) w = 1;
                    else w = 0;
                } else if (ascertain.compare ("polyina")==0) {
                    if (freq[0][j]>0 && freq[0][j] < 1 ) w = 1;
                    else w = 0;
                } else if (ascertain.compare ("hetina")==0) { 
                    a0 --; // Account for ascertainment
                    n0 -= 2;
                    f0 = ((double)a0)/n0;
                    w = 2*freq[0][j]*(1-freq[0][j]);
                }

                double h0 = 0 ;
                double h1 = 0 ;
                
                if (n0 > 1)
                    h0 = a0*(n0-a0)/(n0*(n0-1.));
                if (n1 > 1)
                    h1 = a1*(n1-a1)/(n1*(n1-1.));

                double fstn =  pow((f0-f1),2) - h0/n0 -  h1/n1;
                double fstd  = fstn + h0 + h1;

                fstnum += fstn * w;
                fstdenom += fstd * w;
                chrfstnum[rep-1] += fstn * w;
                chrfstdenom[rep-1] += fstd * w;
                chrw[rep-1] +=  w;
                totalw += w;
            }
        	
            for (int i  = 0 , t = 0 ; i < nsamples; i++)  {
                pi[i] = nucdiv (samples[i], nsnp, snpmat +t );
                h[i] = hfay (samples[i], nsnp, snpmat +t );
                th[i] = thetah (samples[i], nsnp, snpmat +t );

                td [i][rep-1] = tajd(samples[i], nsegsites[i], pi[i]);
                thetaw[i][rep-1] = nsegsites[i]/(2*a1f(samples[i]));
                t += samples[i];
            }
            for (int i = 0 ; i < nsamples; i++)
                delete[] freq[i];
            for (int i = 0 ; i < totalnind; i++)
                delete[] snpmat[i];

			if (printsnplist) { 
				for (int j =0  ; j < nsnp; j++){
					snpfs << snpids[j] << "\t" << chrindex << "\t" << genmap[j] << "\t" << physmap[j] << endl;
				}
			}


		} // One set of MS sims

		if (io::debug >= 2) { 
			cout << "Next MS sim " <<endl;
		}
	}
	snpfs.close ();
	if(inp.rdstate() == ios::eofbit)
		cout << "End of file!\n";
	if(inp.rdstate() == ios::badbit)
		cout << "Fatal I/O error!\n";
	if(inp.rdstate() == ios::failbit)
		cout << "Non-fatal I/O error!\n";
	if(inp.rdstate() == ios::goodbit)
		cout << "No errors!\n";
	inp.close ();	
    divfs.close ();

    double fst = fstdenom>0? fstnum/fstdenom:0;
    double dstat = (abba+baba>0)?(baba-abba)/(baba+abba):0;
    double dstat2 = (abba2+baba2>0)?(baba2-abba2)/(baba2+abba2):0;
    if (l2denom > 0 ) {
        l2en /= l2denom;
        l2yn /= l2denom;
    }
    vector<double> thetawmean (nsamples,0);
    vector<double> thetawsd (nsamples,0);
    vector<double> tajdmean (nsamples,0);
    vector<double> tajdsd (nsamples,0);
    for (int i = 0 ; i < nsamples; i++) { 
        pair<double, double> p = functions::meansd (td[i], reps);
        tajdmean[i] = p.first;
        tajdsd[i] = p.second;

        p = functions::meansd (thetaw[i], reps);
        thetawmean[i] = p.first;
        thetawsd[i] = p.second;

    }

    double divr = divan/divad;
    cout << "dstat: " << dstat << "\t" << "fst: " << fst << "\t";
    cout << "tajd0: " << tajdmean[0] << "\t" << "tajd1: " << tajdmean[1] << "\t";
    cout << "thetaw0: " << thetawmean[0] << "\t" << "thetaw1: " << thetawmean[1] << "\t";
    cout << "dstat2: " << dstat2  ;
    cout << "\tl2en: " << l2en << "\tl2yn: " << l2yn;
    cout << "\tdivr: "  << divr ;
    cout << endl;
    cout << "abba: " << abba << "\t" << "baba: " << baba;
    cout << endl;
	outfs.close ();

    int total1 = 0 ; int total2 = 0;
    for (int i = 0 ; i <= nind[0]; i++) {
        total1 += pop1nd10[i];
        total2 += pop1nd01[i];
    }
    
    jafs.print ("joint.afs");

    ofstream afsout ("ooa.afs") ;
    afsout << "# Count\tTotal\tFrequency\tnd00\tnd01\tnd10\tnd11\tn1\td1\tD-statistic\tS-statistic\tn1/total\td1/total\ttotal"<<endl;
    for (int i =  0 ; i <= nind[0] ; i++) { 
        double j = ((double)i)/nind[0];
        int n = pop1nd10[i]  + pop1nd11[i];
        int d = pop1nd01[i]  + pop1nd11[i];
        int total  = pop1nd00[i] + pop1nd01[i] + pop1nd10[i] + pop1nd11[i];
        double nf = (1.*n)/total;
        double df = (1.*d)/total;
        int s  = pop1nd10[i] - pop1nd01[i];
        double dstat = (1.*s) /(pop1nd10[i]+pop1nd01[i]);
        afsout << i << "\t" << nind[0] << "\t" << j << "\t" << pop1nd00[i] << "\t" << pop1nd01[i] << "\t" << pop1nd10[i] <<"\t" << pop1nd11[i] ;
        afsout << "\t" << n <<"\t" <<d<<"\t"<<dstat <<"\t" <<s <<"\t"<<nf<<"\t"<<df<<"\t" << total<<endl;
    }
    afsout.close ();

    afsout.open ("afr.afs") ;
    afsout << "# Count\tTotal\tFrequency\tnd00\tnd01\tnd10\tnd11\tn1\td1\tD-statistic\tS-statistic\tn1/total\td1/total\ttotal"<<endl;
    for (int i =  0 ; i <= nind[1] ; i++) { 
        double j = ((double)i)/nind[1];
        int n = pop2nd10[i]  + pop2nd11[i];
        int d = pop2nd01[i]  + pop2nd11[i];
        int total  = pop2nd00[i] + pop2nd01[i] + pop2nd10[i] + pop2nd11[i];
        double nf = (1.*n)/total;
        double df = (1.*d)/total;
        int s  = pop2nd10[i] - pop2nd01[i];
        double dstat = (1.*s) /(pop2nd10[i]+pop2nd01[i]);
        afsout << i << "\t" << nind[1] << "\t"<< j << "\t" << pop2nd00[i] << "\t" << pop2nd01[i] << "\t" << pop2nd10[i] <<"\t" << pop2nd11[i] ;
        afsout << "\t" << n <<"\t" <<d<<"\t"<<dstat <<"\t" <<s <<"\t"<<nf<<"\t"<<df<<"\t" << total<<endl;
    }
    afsout.close ();


    vector < vector< double> > data;
    int nstat;
    for (int k = 0 ; k < replicates; k ++) {
        jdenom[k] = totalnsnps - jdenom[k];
        for (int j = 0 ; j < pop1nd00.size(); j++) {
            jpop1nd00[k][j] = pop1nd00[j] - jpop1nd00[k][j];
            jpop1nd01[k][j] = pop1nd01[j] - jpop1nd01[k][j];
            jpop1nd10[k][j] = pop1nd10[j] - jpop1nd10[k][j];
            jpop1nd11[k][j] = pop1nd11[j] - jpop1nd11[k][j];
        }    

        for (int j = 0 ; j < pop2nd00.size(); j++) {
            jpop2nd00[k][j] = pop2nd00[j] - jpop2nd00[k][j];
            jpop2nd01[k][j] = pop2nd01[j] - jpop2nd01[k][j];
            jpop2nd10[k][j] = pop2nd10[j] - jpop2nd10[k][j];
            jpop2nd11[k][j] = pop2nd11[j] - jpop2nd11[k][j];
        }    
        jdivan[k] = divan - jdivan[k];
        jdivad[k] = divad - jdivad[k];
        jdivr[k] = jdivan[k]/jdivad[k];

        ostringstream oss;
        if (io::debug >= 2)
             oss << "jack." << k << ".afs";
        else
            oss << "jack.afs";

        ofstream tmpfs (oss.str());
        tmpfs << "# Count\tTotal\tFrequency\tnd00\tnd01\tnd10\tnd11\tn1\td1\tD-statistic\tS-statistic\tn1/total\td1/total\ttotal"<<endl;
        for (int i =  0 ; i <= nind[1] ; i++) { 
            double j = ((double)i)/nind[1];
            int n = jpop2nd10[k][i]  + jpop2nd11[k][i];
            int d = jpop2nd01[k][i]  + jpop2nd11[k][i];
            int total  = jpop2nd00[k][i] + jpop2nd01[k][i] + jpop2nd10[k][i] + jpop2nd11[k][i];
            double nf = (1.*n)/total;
            double df = (1.*d)/total;
            int s  = jpop2nd10[k][i] - jpop2nd01[k][i];
            double dstat = (1.*s) /(jpop2nd10[k][i]+jpop2nd01[k][i]);
            tmpfs << i << "\t" << nind[1] << "\t"<< j << "\t" << jpop2nd00[k][i] << "\t" << jpop2nd01[k][i] << "\t" << jpop2nd10[k][i] <<"\t" << jpop2nd11[k][i] ;
            tmpfs << "\t" << n <<"\t" <<d<<"\t"<<dstat <<"\t" <<s <<"\t"<<nf<<"\t"<<df<<"\t" << total<<endl;
        }
        tmpfs.close ();
        ostringstream oss1 ;
        if (io::debug >=2  )
            oss1 << "/home/ss545/neander/rolloff/scripts/analyze jack." << k << ".afs 24";
        else
            oss1 << "/home/ss545/neander/rolloff/scripts/analyze jack.afs 24";
        system (oss1.str().c_str());

        if (io::debug >= 2) {
            ostringstream tmposs; tmposs << "cp stats stats."<<k;
            system (tmposs.str().c_str());
        }

        ifstream tmp1fs ("stats");
        string line;
        std::getline ( tmp1fs, line);
        if (line[0] =='#'){
            std::getline ( tmp1fs, line);
        }
        vector<string> toks;
        functions::tokenize (line.c_str(), toks ," \t");
        if (k==0) {
            nstat = toks.size();
            for (int j = 0 ; j < toks.size(); j++) 
                data.push_back (vector<double> ());
        }

        for (int j = 0 ; j < data.size(); j++){  
            double x  = atof(toks[j].c_str());
            data[j].push_back (x);
        }

    }

    if (io::debug >=2 )  {
        cout << "data size = " << data.size()  << endl;
        for (int k  = 0 ; k < replicates; k++){ 
            for (int j = 0 ; j < data.size();j++){
                cout << data[j][k] << "\t";
            }
            cout << endl;
        }
    }

    system ("/home/ss545/neander/rolloff/scripts/analyze afr.afs 24");
    ifstream tmp1fs ("stats");
    std::getline ( tmp1fs, line);
    if (line[0] =='#'){
        std::getline ( tmp1fs, line);
    }
    vector<string> toks;
    functions::tokenize (line.c_str(), toks ," \t");
    vector<double> estimates;
    for (int i =  0 ; i < nstat; i++)
        estimates.push_back (atof(toks[i].c_str()));


    ofstream statfs ("stats.jack");
    for (int i = 0 ; i < data.size(); i++) {
        pair<double, double> p = functions::weightedjack (data[i], jdenom, estimates[i]);
        statfs << estimates[i] <<"\t" << p.first <<"\t" << p.second <<"\t";
    }
    pair<double, double> p = functions::weightedjack ( jdivr, jdenom, divr);
    statfs << divr << "\t" << p.first << "\t" << p.second;

    statfs<<endl;
    statfs.close ();

    ofstream wfs ("weights");
    for (int i = 0 ; i < jdenom.size (); i++)
        wfs << jdenom[i] << endl;
    wfs.close ();

    if (io::debug>=2) { 
        ofstream testfs ("test");
        for (int j = 0 ; j < replicates; j++) { 
            testfs << data[1][j] << "\t" << jdenom[j] << endl;
        }
        testfs.close ();
    }


    ofstream ofs("dstats.jack");
    double dstatnum = (baba-abba); double dstatdenom = baba+abba;

    ofs << "0" << "\t" << dstatnum << "\t" << dstatdenom << endl;
    for (int i = 0 ; i  < reps; i++) {
        chrabba[i] = abba - chrabba[i];
        chrbaba[i] = baba - chrbaba[i];
        chrsnps[i] = totalnsnps - chrsnps[i];
        double dstat = -(chrabba[i]-chrbaba[i])/(chrabba[i]+chrbaba[i]);
        double d1 = chrbaba[i] - chrabba[i];
        double d2 = chrbaba[i] + chrabba[i];
        ofs << i+1 << "\t" << d1 << "\t" << d2 << endl;
//        ofs << dstat << "\t" << chrsnps[i] << endl;
    }
    ofs.close();

    ofs.open("fst.jack");
    for (int i = 0 ; i  < reps; i++) {
        chrfstnum[i] =  fstnum - chrfstnum[i];
        chrfstdenom[i] =  fstdenom - chrfstdenom[i];
        chrw[i] = totalw - chrw[i];
        double chrfst = chrfstnum[i]/chrfstdenom[i];
        ofs <<  chrfst << "\t" << chrw[i] << endl;
    }
    ofs.close();

}

double statsms::nucdiv( int nsam, int segsites, char **list)
{
	int s;
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
   	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	return( pi ) ;
}

/*   thetah - pi   */
double statsms::hfay( int nsam, int segsites, char **list)
{
	int s;
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
   	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ;
		}
	return( -pi ) ;
}

/* Fay's theta_H  */
double statsms::thetah( int nsam, int segsites, char **list)
{
        int s;
        double pi, p1, nd, nnm1  ;

        pi = 0.0 ;

        nd = nsam;
        nnm1 = nd/(nd-1.0) ;
        for( s = 0; s <segsites; s++){
                p1 = frequency('1', s,nsam,list) ;
                pi += p1*p1 ; 
                }
        return( pi*2.0/( nd*(nd-1.0) )  ) ;
}


int statsms::frequency( char allele,int site,int nsam,  char **list)
{
        int i, count=0;
        for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
        return( count);
}


int statsms::segsub( int nsub, int segsites, char **list )
{
	int i, count = 0 , c1 ;

	for(i=0; i < segsites ; i++){
	  c1 = frequency('1',i,nsub, list);
	  if( ( c1 > 0 ) && ( c1 <nsub )  ) count++;
	  }
	return( count ) ;
}



	
int main (int argc, char* argv[]) {
	if (argc > 1) { 
		statsms c = statsms(string(argv[1]));
/*		if (c.mshot)  
    		c.readmshot (c.reps);
		else 
    		c.readms (c.msfile);
            */
        c.readms (c.msfile);
	} else {
		cerr << "Usage: " <<argv[0] << " <config file> "<<endl;
	}
}
