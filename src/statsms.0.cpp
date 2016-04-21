#include "data.h"
#include "io.h"
#include "functions.h"
#include "stringfn.h"
#include "mathfn.h"
#include "statsms.h"
#include "tajd.h"

ofstream outfs;
ofstream disfs;

statsms::statsms (string configfile) { 
	string s;
	d = new data (configfile);
	d->get_int ("debug", io::debug, 0);
	d->get_string ("param", paramfile, true);
	d->get_string ("ms",msfile,true);
	param = new data (paramfile);
	param->get_int ("length", length, 1000000, true);
	param->get_double ("rho", rho, 1.3e-8, true);
	bgrho = rho/(1+(pow(10,2.5)-11)/(log(10)*30));
	param->get_int ("reps", reps, 100, true);
	param->get_string ("samples", s, true);
	io::println ("samples=  " + s ,2);
	vector<string> toks; 
	functions::tokenize (s,  toks, ",");
	for (int i = 0 ; i < toks.size(); i++)
		samples.push_back (atoi (toks[i].c_str()));
	

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
}


/*
void statsms::readmshot (int reps){  
	vector<double> ld0 (nbins) ;
	vector<double> gld0 (nbins) ;
	vector<double> w0 (nbins) ;
	vector<int> lddenom0 (nbins);
	vector<int> wlddenom0 (nbins);
	vector<double> wld0 (nbins) ;
	vector<double> wgld0 (nbins) ;


	vector<double> ld1 (nbins) ;
	vector<double> gld1 (nbins) ;
	vector<double> w1 (nbins) ;
	vector<int> lddenom1 (nbins);
	vector<int> wlddenom1 (nbins);
	vector<double> wld1 (nbins) ;
	vector<double> wgld1 (nbins) ;


	vector<double> ld2 (nbins) ;
	vector<double> gld2 (nbins) ;
	vector<double> w2 (nbins) ;
	vector<int> lddenom2 (nbins);
	vector<int> wlddenom2 (nbins);
	vector<double> wld2 (nbins) ;
	vector<double> wgld2 (nbins) ;

	vector < vector<double> > jgld2(nbins);
	vector < vector<int> > jld2(nbins);

	vector<double> ld3(nbins) ;
	vector<double> gld3 (nbins) ;
	vector<double> w3 (nbins) ;
	vector<int> lddenom3 (nbins);
	vector<int> wlddenom3 (nbins);
	vector<double> wld3 (nbins) ;
	vector<double> wgld3 (nbins) ;

	vector<double> ld23(nbins) ;
	vector<double> gld23 (nbins) ;
	vector<double> w23 (nbins) ;
	vector<int> lddenom23 (nbins);
	vector<int> wlddenom23 (nbins);
	vector<double> wld23 (nbins) ;
	vector<double> wgld23 (nbins) ;

	vector<double> ld4(nbins) ;
	vector<double> gld4 (nbins) ;
	vector<double> w4 (nbins) ;
	vector<int> lddenom4 (nbins);
	vector<int> wlddenom4 (nbins);
	vector<double> wld4 (nbins) ;
	vector<double> wgld4 (nbins) ;

	vector<double> ld34(nbins) ;
	vector<double> gld34 (nbins) ;
	vector<double> w34 (nbins) ;
	vector<int> lddenom34 (nbins);
	vector<int> wlddenom34 (nbins);
	vector<double> wld34 (nbins) ;
	vector<double> wgld34 (nbins) ;
	vector<double> pd (nbins);
	vector<double> pddenom (nbins);

    vector<double> ld5(nbins) ;
	vector<double> gld5 (nbins) ;
	vector<double> w5 (nbins) ;
	vector<int> lddenom5 (nbins);
	vector<int> wlddenom5 (nbins);
	vector<double> wld5 (nbins) ;
	vector<double> wgld5 (nbins) ;



	int nsnp0 = 0 ;
	int nsnp = 0;
	int nsnp1 = 0 ;
	int nsnp2 = 0 ;
	int nsnp3 = 0 ;
	int nsnp4 = 0 ;
	int nsnp5 = 0 ;
	int nsnp23 = 0 ;
	int nsnp34 = 0 ;
	int totalsnps =  0;


	int rep = 0 ;
	int index = 0;
	string id;
	ofstream mapfs ("maperror.txt");
	ofstream snpfs (snpfile.c_str());



	for (int i2 = start ; i2  < start + reps; i2++){
	
	string filename = msfile + "." + tostring(i2);	
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;

	int hsnum;
	vector<int> hsloc;
	vector<double> hsint;
	vector<double> genpos;
	hsloc.push_back(0);
	genpos.push_back(0);
	double gendis = 0;
	double mapl;
	double ermapl;

	while ( std::getline (inp, line)){
		linenum ++;
		io::println ("line = " + line, 2);
		char c = line[0];
		if (c=='#')
			continue;
		
		if (line.empty())
			continue;

		if (line.find ("mshot")==0) {
			string l1 = line;
			int p1 = l1.find ("-v");
			l1 = l1.substr (p1+3,l1.length()-p1-3);
			vector<string> tokens;
			functions::tokenize(l1,tokens," " );
			hsnum = atoi(tokens[0].c_str());
			for (int i3 = 0 ; i3 < hsnum; i3++){
				int l = atoi(tokens[3*i3+1].c_str());
				int r =atoi(tokens[3*i3+2].c_str());
				double inten = atof(tokens[3*i3+3].c_str());
				if( i3==0 ) {
					gendis = l * bgrho;
					genpos.push_back (gendis);
					gendis += (r-l) * inten * bgrho;
					genpos.push_back (gendis);
				} else { 
					double prev = hsloc.back();
					gendis += (l-prev) * bgrho;
					genpos.push_back (gendis);
					gendis += (r-l) *inten*bgrho;
					genpos.push_back (gendis);
				}
				hsloc.push_back (l);
				hsloc.push_back (r);
				hsint.push_back (bgrho);
				hsint.push_back (bgrho*inten);
				if (io::debug >= 2) { 
					cout << "gendis = " <<  l << "\t" << r << "\t" << inten << "\t" <<  gendis << endl;
				}
						
			}
			if (hsloc.back() < length) {
				double prev= hsloc.back();
				hsint.push_back(bgrho);
				hsloc.push_back(length);
				gendis += (length-prev)*bgrho;
				genpos.push_back(gendis);
			}
			mapl = genpos.back();
			ermapl = gsl_ran_gamma (r, mapl * alpha, 1/alpha); 
			cout << "Map length = " << mapl << "\t" << ermapl << endl;

			if (io::debug>=1) {
				cout << "Hotspots" << endl;
				for (int j1 = 0 ;  j1 < hsloc.size(); j1++){
					cout << hsloc[j1] <<"\t" << hsint[j1] << "\t" << genpos[j1] << endl;
				}
			}
		}

		if ( line.find ("segsites:")==0) { 
			string l1 = line.substr(9,line.length()-9);
			nsnp = atoi (l1.c_str());
			io::println ("SNPs = " + tostring(nsnp),2);
		}

		if ( line.find ("positions:")==0) { 
			vector<double> genmap ;
			vector<double> truegenmap ;
			vector<double> physmap;
			vector<string> snpids;
			vector<vector<int> > geno;
			vector<vector<int> > haplo;
			
			string l1 = line.substr(10,line.length()-10);
			//cout << "l1 = " << l1 << endl;
			//
			io::println  ("rep = " + tostring(rep),2);
			//cout << " rep = " << rep  << endl;
			rep++;
			istringstream ss (line.substr(10,line.length()-10));
			int count = 0 ;
			int cur = 0;
			while (!ss.eof()){
				double pos;
				ss>>pos;
			//	cout << pos << endl;
				double ppos = pos * length;
				double gpos = ppos * rho;
				gpos = pos * mapl;
				while (cur < hsloc.size() && hsloc[cur] <=ppos){
					cur++;
				}
				double truegpos = genpos[cur-1] + ((ppos - hsloc[cur-1])/(hsloc[cur]-hsloc[cur-1]))*(genpos[cur]-genpos[cur-1]);
				if (usegenetic==1) { 
					gpos = truegpos;

					if (io::debug>=1) {
						cout << "cur = " << cur << "\t hsloc = " << hsloc[cur] << "\t ppos = " << ppos <<"\t" << gpos << endl;
					}
				} else if (usegenetic == 2) { 
					gpos = pos * ermapl;
				}


				if (count >= nsnp ) 
					break;
				count++;

				genmap.push_back (gpos);
				truegenmap.push_back (truegpos);
				physmap.push_back (ppos);
				id = "snp" + tostring(index);
				snpids.push_back (id);
				index++;
			}
//			int nsnp = physmap.size ();
			if (io::debug >= 1) {
				cout << "Read in " << nsnp << " snps " <<endl;
				ofstream snpfs ("snps.txt");
				for (int i2 = 0 ; i2< nsnp; i2++){
					snpfs <<  snpids[i2] << "\t" << genmap[i2] << "\t" << physmap[i2] << endl;
				}
				snpfs.close ();
			}
			vector<double> freq1 (nsnp,0);
			vector<double> freq2 (nsnp,0);
			vector<double> freq3 (nsnp,0);
			vector<double> freq4 (nsnp,0);
			vector<int> flag0 (nsnp,1);
			vector<int> flag1 (nsnp,1);
			vector<int> flag2 (nsnp,1);
			vector<int> flag3 (nsnp,1);
			vector<int> flag4 (nsnp,1);
			vector<int> flag5 (nsnp,1);
			vector<int> initflag (nsnp,1);

			vector<double> w(nsnp);

			for (int i = 0 ; i < samples[0]/2;i++){
				geno.push_back (vector<int>(nsnp));
				haplo.push_back (vector<int>(nsnp));
				haplo.push_back (vector<int>(nsnp));
			}

			for (int i = 0 ; i < samples[0]/2;i++){
				string s1, s2;
				std::getline (inp, s1);
				std::getline (inp, s2);
				for (int j = 0 ; j < nsnp ; j++) {
					int snp1 = s1[j] - 48 ;
					int snp2 = s2[j] - 48 ;
					haplo[2*i][j] = snp1;
					haplo[1+2*i][j] = snp2;
					geno[i][j] = snp1+snp2;
					
					freq1[j] += geno[i][j];
				}	
			}
			if (io::debug>=2) {
				for (int i = 0 ; i < samples[0]/2;i++){
					for (int j = 0 ; j < nsnp ; j++) {
						cout << i<<","<<j<<"\t"<<haplo[2*i][j]<<","<<haplo[2*i+1][j]<<","<<geno[i][j]<<endl;
						if (geno[i][j] < 0 || geno[i][j] > 2){ 
							cout << "Error at individual " << i <<",snp " <<j << ", geno = " << geno[i][j]<<endl;
							exit(1);
						}
					}
				}
			}

			for (int j  = 0 ; j < nsnp; j++) 
				freq1[j]/=samples[0];


			for (int i  = 0 ; i < samples[1] ; i++) {
				string s;
				std::getline (inp, s);
				for (int j = 0 ; j < nsnp ; j++) {
					int snp = s[j] - 48;
					freq2[j] += snp;
				}
			}
			for (int j  = 0 ; j < nsnp; j++) 
				freq2[j]/=samples[1];

			for (int i  = 0 ; i < samples[2] ; i++) {
				string s;
				std::getline (inp, s);
				for (int j = 0 ; j < nsnp ; j++) {
					int snp = s[j] - 48;
					freq3[j] += snp;
				}
			}
			for (int j  = 0 ; j < nsnp; j++) 
				freq3[j]/=samples[2];

			if (nsamples == 4) {
				for (int i  = 0 ; i < samples[3] ; i++) {
					string s;
					std::getline (inp, s);
					for (int j = 0 ; j < nsnp ; j++) {
						int snp = s[j] - 48;
						freq4[j] += snp;
					}
				}
				for (int j  = 0 ; j < nsnp; j++) 
					freq4[j]/=samples[3];

			}

			for (int j  = 0 ; j < nsnp; j++)  {
				if (uniform)
					w[j] = 1;
				else 
					w[j] = (freq2[j] - freq3[j]);

				int derived1 = freq2[j]>0;
				int derived2 = freq3[j]>0;

				if (kg) {
					double callp[] = {0.25, 0.5, 0.75, 0.80, 0.9,0.95,0.96,0.97,0.98,0.99}; 
					int n1 = mathfn::round(freq1[j]*samples[0]);
					int v1 = samples[0] - n1;
					int n2 = mathfn::round(freq2[j]*samples[1]);
					int v2 = samples[1] - n2;
					
					double x = gsl_ran_flat ( r, 0,1);
					double t = 1;
					if (n1>=10 && v1 >=10){ 
						t = 0.99;
					} else if (n1==0 || v1==0){
						// Don't do anything here. Depends on whether this is polymorphic.
					} else {
						int ind = n1<v1?n1:v1;
						t = callp[ind-1];
					}
					
					if (n1==0 || v1==0) { 
						initflag[j] = 0;
					} else { 
						if (x > t) 
							flag0[j] = flag1[j] = flag2[j] = flag3[j] = flag4[j] = flag5[j] = initflag[j] = 0;
						else 
							totalsnps ++;
					}


					
				} else { 
					if (freq1[j]>0)
						totalsnps ++;
				}

				switch (ascertain) { 
				case 0:
					break;
				case 1:
					break;
				case 2:
					if (nsamples == 3) { 
                        // derived1 -- derived in Y ?
                        // derived2 -- derived in N ?
                        // N derived, YRI all ancestral
						if ( derived2 == 1 && derived1 == 0 ) {}
						else flag3[j] = 0;
	
                        // N derived
						if (derived2==1){}
						else flag1[j] = 0;


                        // N derived, CEU < cutoff
						if (freq1[j] > fcutoff)
							flag2[j] = 0;
						if (derived2==0)
							flag2[j] = 0;

                        // N derived, YRI all ancestral
                        // + symmetric
						if (derived2==1){
							if ((freq1[j]>0 && derived1==0) || (freq1[j]==0 && freq2[j]>0)){}
							else
								flag4[j] = 0;
//							if ((freq1[j]>0 && freq2[j]==0) || (freq1[j]==0 && freq2[j]==0)){}
//							else
//								flag4[j] = 0;
						} else {
							flag4[j] = 0;
						}



						if(io::debug >= 2)
							cout << "j= " <<j <<"\tflag = " << flag2[j]<<endl;

					} else {
						if ( asn ) {
							int derived3 = freq4[j] > 0;
							if ( derived2 == 1 && derived1 == 0 ) {}
							else flag3[j] = 0;

							if (derived2==1){}
							else flag1[j] = 0;


							if (freq1[j] > fcutoff)
								flag2[j] = 0;
							if (derived2==0)
								flag2[j] = 0;


						} else { 
							int derived3 = freq4[j];

							if ( derived2 == 1 && derived3==0 && derived1 == 0 ) {}
							else flag3[j] = 0;

							if (derived2==1 && derived3==0){}
							else flag1[j] = 0;

							if (freq1[j] > fcutoff)
								flag2[j] = 0;
							if (derived2==1 && derived3==0) {}
							else
								flag2[j] = 0;

                            if (freq1[j] > fcutoff)
                                flag5[j] = 0;
                            if (derived2==0 && derived3==1){}
                            else 
                                flag5[j] = 0;

							if(io::debug >= 2)
								cout << "j= " <<j <<"\tflag = " << flag2[j]<<endl;


							if (derived2==1){
								if ((freq1[j]>0 && freq2[j]==0) || (freq1[j]==0 && freq2[j]>0)){}
								else
									flag4[j] = 0;
							} else {
								flag4[j] = 0;
							}

						}

					}

				}

				if (asn) {
					if (poly==1) {
						if (freq1[j]==0 && freq4[j]==0){
							flag2[j] = flag3[j] =  0;
							flag0[j] = flag1[j] =  0;
						} 
					}
				} else {
					if (poly==1 && (freq1[j]==0) ){
						flag2[j] = flag3[j] =  0;
						flag0[j] = flag1[j] =  0;
                        flag5[j] = 0;
					} 
				}

				if (flag0[j])
					nsnp0++;
				if (flag1[j])
					nsnp1++;
				if (flag2[j])
					nsnp2++;
				if (flag5[j])
					nsnp5++;
				if (flag3[j])
					nsnp3++;
				if (flag2[j] && flag3[j])
					nsnp23++;


				if (flag2[j] && io::debug>=1) { 
					cout << "Accepted =  " << snpids[j] << "\t" << freq1[j] << "\t" << freq2[j] << "\t" << freq3[j]<<endl;
				}
			}

			if (io::debug >= 2) {
				int c1 = 0 ; int c2 = 0; int  c3 = 0; int c0 = 0 ;
                int c5 = 0;
				cout << "nsnp = " << nsnp << endl;
				for (int j = 0 ; j < nsnp ; j++){
					c1 += flag1[j];
					c0 += flag0[j];
					c2 += flag2[j];
					c3 += flag3[j];
                    c5 += flag5[j];
				}
				cout << c0 << "," << c1 << "," << c2 << "," << c3 << "," << endl;
			}
			
			if (usegenetic > 0 ) { 
				double prev ;
				double prevpd;
				bool flag =  false;
				for (int j = 0; j < nsnp ; j++){ 
					if (!flag) {
						if (initflag[j]){
							flag = true;
							prev = genmap[j];
							prevpd = physmap[j];
							double badpos = (prevpd/length)*mapl;
							snpfs << "chr"<<i2<<":"<<j<<"\t"<<i2 << "\t" <<  genmap[j] << "\t" << truegenmap[j] << endl;
						}
					} else { 
						if (initflag[j]){
							double dis = genmap[j] - prev;
							double pdis = physmap[j] - prevpd;
							double obdis = (pdis/length)*mapl;
							prev = genmap[j];
							prevpd = physmap[j];
							mapfs << i2 << "\t" << obdis << "\t" << dis << endl;
							double badpos = (physmap[j]/length)*mapl;
							snpfs << "chr"<<i2<<":"<<j<<"\t"<<i2 << "\t" <<  genmap[j] << "\t" << truegenmap[j] << endl;

						}
					}	
				}
			} else { 
				for (int j = 0; j < nsnp ; j++){ 
					snpfs << "chr"<<i2<<":"<<j<<"\t"<<i2 << "\t" <<  genmap[j] << "\t" << truegenmap[j] << endl;
				}
			}

			for (int j =  0 ; j < nsnp; j++) {
				for (int k = j+1; k < nsnp ; k++){ 
					double dis = genmap[k]-genmap[j];
					double pdis = physmap[k]  - physmap[j];
					if ( dis >= maxdis)
						continue;
					
						
//					if ((flag2[j]==1 && flag2[k]==1)||(flag3[j]==1 && flag3[k]==1) ||(flag1[j]==1&&flag1[k]==1) || (flag0[j]==1 && flag0[k]==1)) {}
//					else
//						continue;
						
					
					if (fast)  {	
						if ((flag2[j]==1 && flag2[k]==1)||(flag3[j]==1 && flag3[k]==1)||(flag4[j]&&flag4[k]) || (flag5[j]==1 && flag5[k]==1)) {}
						else
							continue;
					} else { 
						if ((flag2[j]==1 && flag2[k]==1)||(flag3[j]==1 && flag3[k]==1) ||(flag1[j]==1&&flag1[k]==1) || (flag0[j]==1 && flag0[k]==1)|| (flag5[j]==1 && flag5[k]==1)) {}
						else
							continue;
					}

					int ind =(int)( dis/binsize);
					pd[ind] += pdis;
					pddenom[ind]++;

					if (io::debug >= 3) {
						cout << "samples = " << samples[0] <<"," << samples[1]<<","<<samples[2]<<endl;
					}
					// Haplotype based
					double a[2][2];
					a[0][0] = a[0][1] = a[1][0] = a[1][1]= 0 ;
					for (int i  = 0 ;i < samples[0]; i++)
						a[haplo[i][j]][haplo[i][k]]++;
					double d1 = (a[1][1]*a[0][0] - a[0][1] * a[1][0])/(samples[0]*samples[0]);

					// genotype based
					double b[3][3];
					for (int i1  = 0; i1 < 3; i1++) {
						for (int j1 = 0 ; j1 < 3; j1++){
							b[i1][j1] = 0 ;
						}
					}


					for (int i  = 0 ;i < samples[0]/2; i++)
						b[geno[i][j]][geno[i][k]]++;

//					continue;

					double d2 = 0 ;
					double p  = 0;
					double q = 0;
					for (int i1  = 0; i1 < 3; i1++){
						for (int j1 = 0 ; j1 < 3; j1++){
							d2 += i1*j1*b[i1][j1];
							p += i1*b[i1][j1];
							q += j1*b[i1][j1];
						}
					}	
					d2 /= (samples[0]/2);
					d2-=4*p*q/(samples[0]*samples[0]);

//					continue;
                    if (flag0[j] && flag0[k]) {
                        ld0[ind] += d1;
                        gld0[ind] += d2;
                        lddenom0[ind]++;

                        double wj0, wk0;
                        if (freq3[j]-freq2[j] > 0.9)
                            wj0 = freq3[j] - freq2[j];
                        else if (freq3[j]-freq2[j] < 0.9)
                            wj0 = -(freq3[j] - freq2[j]);
                        else 
                            wj0 = 0;

                        if (freq3[k]-freq2[k] > 0.9)
                            wk0 = freq3[k] - freq2[k];
                        else if (freq3[k]-freq2[k] < 0.9)
                            wk0 = -(freq3[k] - freq2[k]);
                        else 
                            wk0 = 0;


                        wld0[ind] += wj0*wk0*d1;
                        wgld0[ind] += wj0*wk0*d2;
                        wlddenom0[ind] += (wj0*wk0)*(wj0*wk0);

                    }



                    if (flag1[j] && flag1[k]) {
                        ld1[ind] += d1;
                        gld1[ind] += d2;
                        lddenom1[ind]++;

                        wld1[ind] += w[j]*w[k]*d1;
                        wgld1[ind] += w[j]*w[k]*d2;
                        wlddenom1[ind] += (w[j]*w[k])*(w[j]*w[k]);

                    }

                    if (flag2[j] && flag2[k]){
                        ld2[ind] += d1;
                        gld2[ind] += d2;
                        lddenom2[ind]++;

                        wld2[ind] += w[j]*w[k]*d1;
                        wgld2[ind] += w[j]*w[k]*d2;
                        wlddenom2[ind] += (w[j]*w[k])*(w[j]*w[k]);
                    }

                    if (flag5[j] && flag5[k]){
                        ld5[ind] += d1;
                        gld5[ind] += d2;
                        lddenom5[ind]++;

                        wld5[ind] += w[j]*w[k]*d1;
                        wgld5[ind] += w[j]*w[k]*d2;
                        wlddenom5[ind] += (w[j]*w[k])*(w[j]*w[k]);
                    }

                    if (flag3[j] && flag3[k]){
                        ld3[ind] += d1;
                        gld3[ind] += d2;
                        lddenom3[ind]++;

                        wld3[ind] += w[j]*w[k]*d1;
                        wgld3[ind] += w[j]*w[k]*d2;
                        wlddenom3[ind] += (w[j]*w[k])*(w[j]*w[k]);

                        if (flag2[j] && flag2[k]) { 
                            ld23[ind] += d1;
                            gld23[ind] += d2;
                            lddenom23[ind]++;

                            wld23[ind] += w[j]*w[k]*d1;
                            wgld23[ind] += w[j]*w[k]*d2;
                            wlddenom23[ind] += (w[j]*w[k])*(w[j]*w[k]);
                        }
                    }

                    if (flag4[j] && flag4[k]) {
                        ld4[ind] += d1;
                        gld4[ind] += d2;
                        lddenom4[ind]++;

                        wld4[ind] += w[j]*w[k]*d1;
                        wgld4[ind] += w[j]*w[k]*d2;
                        wlddenom4[ind] += (w[j]*w[k])*(w[j]*w[k]);

                        if (freq1[j]==0 || freq1[j]<=fcutoff) { 
                            ld34[ind] += d1;
                            gld34[ind] += d2;
                            lddenom34[ind]++;

                            wld34[ind] += w[j]*w[k]*d1;
                            wgld34[ind] += w[j]*w[k]*d2;
                            wlddenom34[ind] += (w[j]*w[k])*(w[j]*w[k]);
                        }

                    }


				}
			}


		} // One set of MS sims

		if (io::debug >= 2) { 
			cout << "Next MS sim " <<endl;
		}

	}

	if(inp.rdstate() == ios::eofbit)
		cout << "End of file!\n";
	if(inp.rdstate() == ios::badbit)
		cout << "Fatal I/O error!\n";
	if(inp.rdstate() == ios::failbit)
		cout << "Non-fatal I/O error!\n";
	if(inp.rdstate() == ios::goodbit)
		cout << "No errors!\n";
	inp.close ();	

	}
	mapfs.close ();
	snpfs.close ();

	for (int i = 0 ; i < nbins; i++) {
		if (pddenom[i]>0) {
			pd[i]/=pddenom[i];
		}

		if (lddenom0[i]>0) { 
			ld0[i]/=lddenom0[i];
			gld0[i]/=lddenom0[i];
		}
		if (wlddenom0[i]>0) { 
			wld0[i]/=wlddenom0[i];
			wgld0[i]/=wlddenom0[i];
		}
		if (lddenom1[i]>0) { 
			ld1[i]/=lddenom1[i];
			gld1[i]/=lddenom1[i];
		}
		if (wlddenom1[i]>0) { 
			wld1[i]/=wlddenom1[i];
			wgld1[i]/=wlddenom1[i];
		}

		if (lddenom2[i]>0) { 
			ld2[i]/=lddenom2[i];
			gld2[i]/=lddenom2[i];
		}
		if (wlddenom2[i]>0) { 
			wld2[i]/=wlddenom2[i];
			wgld2[i]/=wlddenom2[i];
		}

		if (lddenom5[i]>0) { 
			ld5[i]/=lddenom5[i];
			gld5[i]/=lddenom5[i];
		}
		if (wlddenom5[i]>0) { 
			wld5[i]/=wlddenom5[i];
			wgld5[i]/=wlddenom5[i];
		}

		if (lddenom3[i]>0) { 
			ld3[i]/=lddenom3[i];
			gld3[i]/=lddenom3[i];
		}
		if (wlddenom3[i]>0) { 
			wld3[i]/=wlddenom3[i];
			wgld3[i]/=wlddenom3[i];
		}
		if (lddenom23[i]>0) { 
			ld23[i]/=lddenom23[i];
			gld23[i]/=lddenom23[i];
		}
		if (wlddenom23[i]>0) { 
			wld23[i]/=wlddenom23[i];
			wgld23[i]/=wlddenom23[i];
		}
		if (lddenom4[i]>0) { 
			ld4[i]/=lddenom4[i];
			gld4[i]/=lddenom4[i];
		}
		if (wlddenom4[i]>0) { 
			wld4[i]/=wlddenom4[i];
			wgld4[i]/=wlddenom4[i];
		}
		if (lddenom34[i]>0) { 
			ld34[i]/=lddenom34[i];
			gld34[i]/=lddenom34[i];
		}
		if (wlddenom34[i]>0) { 
			wld34[i]/=wlddenom34[i];
			wgld34[i]/=wlddenom34[i];
		}


		double j = (i+1.0)/nbins;
		if (out) { 
            // N derived , freq < cutoff
			outfs << j <<"\t" << ld2[i] << "\t" << gld2[i] << "\t" << wld2[i] << "\t" << wgld2[i] << "\t" << lddenom2[i] ;
            // N derived, YRI all ancestral
			outfs << "\t" << ld3[i] << "\t" << gld3[i] << "\t" << wld3[i] << "\t" << wgld3[i] << "\t" << lddenom3[i] ;
            // Intersect 2 and 3
			outfs << "\t" << ld23[i] << "\t" << gld23[i] << "\t" << wld23[i] << "\t" << wgld23[i] << "\t" << lddenom23[i] ;
            // N derived
			outfs << "\t" << ld1[i] << "\t" << gld1[i] << "\t" << wld1[i] << "\t" << wgld1[i] << "\t" << lddenom1[i] ;
            // polymorphic
			outfs << "\t" << ld0[i] << "\t" << gld0[i] << "\t" << wld0[i] << "\t" << wgld0[i] << "\t" << lddenom0[i] ;
            // N derived, YRI all ancestral
            // + symmetric
			outfs << "\t" << ld4[i] << "\t" << gld4[i] << "\t" << wld4[i] << "\t" << wgld4[i] << "\t" << lddenom4[i] ;
            // N derived, YRI all ancestral, CEU freq <= fcutoff
            // N derived, YRI poly, CEU==0
			outfs << "\t" << ld34[i] << "\t" << gld34[i] << "\t" << wld34[i] << "\t" << wgld34[i] << "\t" << lddenom34[i] ;
            // N ancestral, X derived , freq < cutoff
			outfs << j <<"\t" << ld5[i] << "\t" << gld5[i] << "\t" << wld5[i] << "\t" << wgld5[i] << "\t" << lddenom5[i] ;
			outfs << endl;
			disfs << j << "\t" << pd[i] << "\t" << pddenom[i] << endl;

        } else { 
            cout << j <<"\t" << ld2[i] << "\t" << gld2[i] << "\t" << wld2[i] << "\t" << wgld2[i] << "\t" << lddenom2[i] ;
            cout << "\t" << ld3[i] << "\t" << gld3[i] << "\t" << wld3[i] << "\t" << wgld3[i] << "\t" << lddenom3[i] ;
            cout << "\t" << ld23[i] << "\t" << gld23[i] << "\t" << wld23[i] << "\t" << wgld23[i] << "\t" << lddenom23[i] ;
            cout << "\t" << ld1[i] << "\t" << gld1[i] << "\t" << wld1[i] << "\t" << wgld1[i] << "\t" << lddenom1[i] ;
            cout << "\t" << ld0[i] << "\t" << gld0[i] << "\t" << wld0[i] << "\t" << wgld0[i] << "\t" << lddenom0[i] ;
            cout << "\t" << ld4[i] << "\t" << gld4[i] << "\t" << wld4[i] << "\t" << wgld4[i] << "\t" << lddenom4[i] ;
            cout << "\t" << ld34[i] << "\t" << gld34[i] << "\t" << wld34[i] << "\t" << wgld34[i] << "\t" << lddenom34[i] ;
            // N ancestral, X derived , freq < cutoff
            cout << j <<"\t" << ld5[i] << "\t" << gld5[i] << "\t" << wld5[i] << "\t" << wgld5[i] << "\t" << lddenom5[i] ;
            cout << endl;
        }
	}
	cout << "Total number of snps =  " << totalsnps << endl;
	if (io::debug >= 1) { 
		cout << "Number of snps =  " << nsnp2 << "\t" <<nsnp3 <<"\t"<<nsnp23<<endl;
	}
	outfs.close ();
	disfs.close ();
}*/


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

    vector<int> pop1nd10 (nind[0]+1,0);
    vector<int> pop1nd01 (nind[0]+1,0);
    vector<int> pop2nd10 (nind[1]+1,0);
    vector<int> pop2nd01 (nind[1]+1,0);
    vector<int> pop1nd00 (nind[0]+1,0);
    vector<int> pop2nd00 (nind[1]+1,0);
    vector<int> pop1nd11 (nind[0]+1,0);
    vector<int> pop2nd11 (nind[1]+1,0);
    vector<int> pop1n1 (nind[0]+1,0);
    vector<int> pop2n1 (nind[0]+1,0);
    vector<int> pop1d1 (nind[0]+1,0);
    vector<int> pop2d1 (nind[1]+1,0);
    vector<int> pop1n0 (nind[0]+1,0);
    vector<int> pop2n0 (nind[0]+1,0);
    vector<int> pop1d0 (nind[0]+1,0);
    vector<int> pop2d0 (nind[1]+1,0);

    vector<double> chrabba = vector<double>(reps,0);
    vector<double> chrbaba = vector<double>(reps,0);
    vector<double> chrsnps = vector<double>(reps,0);
    vector<double> chrfstnum = vector<double> (reps,0);
    vector<double> chrfstdenom = vector<double> (reps,0);
    vector<double> chrw = vector<double> (reps,0);
    int totalnsnps = 0;

	int rep = 0 ;
	int index = 0;
	string id;

	while ( std::getline (inp, line)){
		linenum ++;
		io::println ("line = " + line, 2);
		char c = line[0];
		if (c=='#')
			continue;
		
		if (line.empty())
			continue;

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
			io::println  ("rep = " + tostring(rep),2);
			//cout << " rep = " << rep  << endl;
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
            

            chrsnps[rep-1] = nsnp;
            totalnsnps += nsnp;
            for (int j = 0 ; j < nsnp ; j++)  {

                if (nsamples >= 5 ) {
                    if (ascertain.compare ("african")==0) { 
                        if (freq[2][j] == 0  && freq[3][j] == 0 && freq[4][j] ==0 && freq[1][j]==0)
                            continue;
                    }
                }
                if (freq[2][j]>=1){
                    abba += freq[1][j] * ( 1-freq[0][j]);
                    baba += freq[0][j] * ( 1-freq[1][j]);
                    chrabba[rep-1] += freq[1][j] * ( 1-freq[0][j]);
                    chrbaba [rep-1] += freq[0][j] * ( 1-freq[1][j]);
                }
                if (freq[2][j]>0){
                    int ind = round(freq[0][j]*samples[0]);
                    pop1n1[ind]++;
                    ind = round(freq[1][j]*samples[1]);
                    pop2n1[ind]++;
                } else  {
                    int ind = round(freq[0][j]*samples[0]);
                    pop1n0[ind]++;
                    ind = round(freq[1][j]*samples[1]);
                    pop2n0[ind]++;
                }
                if (nsamples >= 4 ) { 
                    if (freq[3][j]>=1){
                        abba2 += freq[1][j] * ( 1-freq[0][j]);
                        baba2 += freq[0][j] * ( 1-freq[1][j]);
                    }
                    if (freq[3][j] >  0 ) {
                        int ind = round(freq[0][j]*samples[0]);
                        pop1d1[ind]++;
                        ind = round(freq[1][j]*samples[1]);
                        pop2d1[ind]++;
                    } else {
                        int ind = round(freq[0][j]*samples[0]);
                        pop1d0[ind]++;
                        ind = round(freq[1][j]*samples[1]);
                        pop2d0[ind]++;
                    }
                   //  cout << "here1" << endl;
                    if (freq[2][j]>0 && freq[3][j]==0) {
                        int ind = round(freq[0][j]*samples[0]);
                        pop1nd10[ind]++;
                    //    cout << "here2\t" << ind  << "\t" << freq[0][j] << "\t" << samples[0]<< endl;
                        ind = round(freq[1][j]*samples[1]);
                        pop2nd10[ind]++;
                    }
                    if (freq[3][j]>0 && freq[2][j]==0) {
                        int ind = round(freq[0][j]*samples[0]);
                        pop1nd01[ind]++;
                    //    cout << "here3\t" << ind  << "\t" << freq[0][j] << "\t" << samples[0]<< endl;
                        ind = round(freq[1][j]*samples[1]);
                        pop2nd01[ind]++;
                    }
                    if (freq[2][j]==0 && freq[3][j]==0) {
                        int ind = round(freq[0][j]*samples[0]);
                        pop1nd00[ind]++;
                        ind = round(freq[1][j]*samples[1]);
                        pop2nd00[ind]++;
                    }
                    if (freq[2][j]>0 && freq[3][j]>0) {
                        int ind = round(freq[0][j]*samples[0]);
                        pop1nd11[ind]++;
                        ind = round(freq[1][j]*samples[1]);
                        pop2nd11[ind]++;
                    }
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

    cout << "dstat: " << dstat << "\t" << "fst: " << fst << "\t";
    cout << "tajd0: " << tajdmean[0] << "\t" << "tajd1: " << tajdmean[1] << "\t";
    cout << "thetaw0: " << thetawmean[0] << "\t" << "thetaw1: " << thetawmean[1] << "\t";
    cout << "dstat2: " << dstat2  ;
    cout << "\tl2en: " << l2en << "\tl2yn: " << l2yn;
    cout << endl;
	outfs.close ();

    int total1 = 0 ; int total2 = 0;
    for (int i = 0 ; i <= nind[0]; i++) {
        total1 += pop1nd10[i];
        total2 += pop1nd01[i];
    }
    ofstream afs ("pop1.afs");
    afs << "# Total number of sites = " << total1 << "\t" << total2 << endl;
    int diff = pop1nd10[0]-pop1nd01[0];
    double d = ((double)diff)/(pop1nd10[0]+pop1nd01[0]);
    afs << "#0\t0.00\t" << pop1nd10[0] << "\t" << pop1nd01[0]  <<"\t"<< diff << "\t"<<d <<endl;
    double j = ((double)nind[0])/nind[0];
    diff = pop1nd10[nind[0]]-pop1nd01[nind[0]];
    d = ((double)diff)/(pop1nd10[nind[0]]+pop1nd01[nind[0]]);
    afs << "#"<<nind[0]<<"\t"<<j<<"\t" << pop1nd10[nind[0]] << "\t" << pop1nd01[nind[0]] <<"\t"<<diff << "\t"<<d << endl;
    for (int i = 1 ; i < nind[0]; i++) {
        double j = ((double)i)/nind[0];
        int diff = pop1nd10[i]-pop1nd01[i];
        double  d = ((double)diff)/ (pop1nd10[i]+pop1nd01[i]);
        afs << i << "\t" << j <<"\t"<<pop1nd10[i] << "\t"<<pop1nd01[i] <<"\t" << diff <<"\t" <<d << endl;
    }
    afs.close ();

    total1 = 0 ; total2 = 0;
    for (int i = 0 ; i <= nind[1]; i++) {
        total1 += pop2nd10[i];
        total2 += pop2nd01[i];
    }
    afs.open ("pop2.afs");
    afs << "# Total number of sites = " << total1 << "\t" << total2 << endl;
    diff = pop2nd10[0]-pop2nd01[0];
    d = ((double)diff)/(pop2nd10[0]+pop2nd01[0]);
    afs << "#0\t0.00\t" << pop2nd10[0] << "\t" << pop2nd01[0]  <<"\t"<< diff << "\t"<<d <<endl;
    j = ((double)nind[0])/nind[0];
    diff = pop2nd10[nind[0]]-pop2nd01[nind[0]];
    d = ((double)diff)/(pop2nd10[nind[0]]+pop2nd01[nind[0]]);
    afs << "#"<<nind[0]<<"\t"<<j<<"\t" << pop2nd10[nind[0]] << "\t" << pop2nd01[nind[0]] <<"\t"<<diff << "\t"<<d << endl;

    for (int i = 1 ; i < nind[1]; i++) {
        double j = ((double)i)/nind[0];
        int diff = pop2nd10[i]-pop2nd01[i];
        double  d = ((double)diff)/ (pop2nd10[i]+pop2nd01[i]);
        afs << i << "\t" << j <<"\t"<<pop2nd10[i] << "\t"<<pop2nd01[i] <<"\t" << diff << "\t" <<d <<endl;
    }
    afs.close ();

    total1 = 0 ; total2 = 0;
    for (int i = 0 ; i <= nind[0]; i++) {
        total1 += pop1n1[i];
    }
    for (int i = 0 ; i <= nind[1]; i++) {
        total2 += pop2n1[i];
    }

    afs.open ("n.afs");
    afs << "# Total number of sites = " << total1 << "\t" << total2 << endl;
    afs << "# Fixed ancestral " << pop1n1[0] << "\t" << pop2n1[0] << "\t" << (1.*pop1n1[0])/(pop1n1[0]+pop1n0[0]) << "\t" << (1.*pop2n1[0])/(pop2n1[0]+pop2n0[0]) << endl;
    afs << "# Fixed derived " << pop1n1[nind[0]] << "\t" << pop2n1[nind[0]] << "\t" << (1.*pop1n1[nind[0]])/(pop1n1[nind[0]]+pop1n0[nind[0]])<<"\t"<<(1.*pop2n1[nind[0]])/(pop2n1[nind[0]]+pop2n0[nind[0]])<< endl;
    for (int i  = 1 ; i < nind[0]; i++){
        double j = ((double)i)/nind[0];
        double x = (1.*pop1n1[i])/(pop1n1[i]+pop1n0[i]);
        double y = (1.*pop2n1[i])/(pop2n1[i]+pop2n0[i]);
        afs << i << "\t" << j <<"\t"<<pop1n1[i] << "\t"<<pop2n1[i] << "\t" << x << "\t" << y << endl;
    }
    afs.close ();

    afs.open ("d.afs");
    afs << "# Total number of sites = " << total1 << "\t" << total2 << endl;
    afs << "# Fixed ancestral " << pop1d1[0] << "\t" << pop2d1[0] << "\t" << (1.*pop1d1[0])/(pop1d1[0]+pop1d0[0]) << "\t" << (1.*pop2d1[0])/(pop2d1[0]+pop2d0[0]) << endl;
    afs << "# Fixed derived " << pop1d1[nind[0]] << "\t" << pop2d1[nind[0]] << "\t" << (1.*pop1d1[nind[0]])/(pop1d1[nind[0]]+pop1d0[nind[0]])<<"\t"<<(1.*pop2d1[nind[0]])/(pop2d1[nind[0]]+pop2d0[nind[0]])<< endl;
    for (int i  = 1 ; i < nind[0]; i++){
        double j = ((double)i)/nind[0];
        double x = (1.*pop1d1[i])/(pop1d1[i]+pop1d0[i]);
        double y = (1.*pop2d1[i])/(pop2d1[i]+pop2d0[i]);
        afs << i << "\t" << j <<"\t"<<pop1d1[i] << "\t"<<pop2d1[i] << "\t" << x << "\t" << y << endl;
    }
    afs.close ();


    afs.open ("pop1.total.afs") ;
    afs << "# Count\tTotal\tFrequency\tnd00\tnd01\tnd10\tnd11\tn1\td1\tD-statistic\tS-statistic\tn1/total\td1/total\ttotal"<<endl;
    for (int i =  0 ; i <= nind[0] ; i++) { 
        double j = ((double)i)/nind[0];
        int n = pop1nd10[i]  + pop1nd11[i];
        int d = pop1nd01[i]  + pop1nd11[i];
        int total  = pop1nd00[i] + pop1nd01[i] + pop1nd10[i] + pop1nd11[i];
        double nf = (1.*n)/total;
        double df = (1.*d)/total;
        int s  = pop1nd10[i] - pop1nd01[i];
        double dstat = (1.*s) /(pop1nd10[i]+pop1nd01[i]);
        afs << i << "\t" << nind[0] << "\t" << j << "\t" << pop1nd00[i] << "\t" << pop1nd01[i] << "\t" << pop1nd10[i] <<"\t" << pop1nd11[i] ;
        afs << "\t" << n <<"\t" <<d<<"\t"<<dstat <<"\t" <<s <<"\t"<<nf<<"\t"<<df<<"\t" << total<<endl;
    }
    afs.close ();

    afs.open ("pop2.total.afs") ;
    afs << "# Count\tTotal\tFrequency\tnd00\tnd01\tnd10\tnd11\tn1\td1\tD-statistic\tS-statistic\tn1/total\td1/total\ttotal"<<endl;
    for (int i =  0 ; i <= nind[1] ; i++) { 
        double j = ((double)i)/nind[1];
        int n = pop2nd10[i]  + pop2nd11[i];
        int d = pop2nd01[i]  + pop2nd11[i];
        int total  = pop2nd00[i] + pop2nd01[i] + pop2nd10[i] + pop2nd11[i];
        double nf = (1.*n)/total;
        double df = (1.*d)/total;
        int s  = pop2nd10[i] - pop2nd01[i];
        double dstat = (1.*s) /(pop2nd10[i]+pop2nd01[i]);
        afs << i << "\t" << nind[1] << "\t"<< j << "\t" << pop2nd00[i] << "\t" << pop2nd01[i] << "\t" << pop2nd10[i] <<"\t" << pop2nd11[i] ;
        afs << "\t" << n <<"\t" <<d<<"\t"<<dstat <<"\t" <<s <<"\t"<<nf<<"\t"<<df<<"\t" << total<<endl;
    }
    afs.close ();



    ofstream ofs("dstats.jack");
    for (int i = 0 ; i  < reps; i++) {
        chrabba[i] = abba - chrabba[i];
        chrbaba[i] = baba - chrbaba[i];
        chrsnps[i] = totalnsnps - chrsnps[i];
        double dstat = -(chrabba[i]-chrbaba[i])/(chrabba[i]+chrbaba[i]);
        ofs << dstat << "\t" << chrsnps[i] << endl;
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
