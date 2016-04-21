#include "data.h"
#include "io.h"
#include "functions.h"
#include "stringfn.h"
#include "mathfn.h"

ofstream outfs;

class twopointsim  {
	public:
		string paramfile;
		string msfile;
		data *param;

		twopointsim (string);
		void run();
		void runmany ();

		int gen;
		int popsize;
		double recrate;
		bool varyrec;
		int bstart;
		int bend;
		int bsize;
		bool bkpt;
		int ensize;
		int seed;

		vector<double> initfreq;


		const gsl_rng_type * T;
	        gsl_rng * r;
};


twopointsim::twopointsim (string paramfile) { 
	param = new data (paramfile);
	param->get_int ("generations",gen,2000,false);
	param->get_int ("popsize",popsize,10000,false);
	param->get_double ("recrate",recrate,0.1,false);
	recrate*=1e-2;
	param->get_int ("bstart", bstart,0,false);
	param->get_int ("bend", bend,0,false);
	param->get_int ("bsize", bsize,0,false);
	param->get_int ("ensize", ensize, 10000, false);
	param->get_int ("seed",seed,1,false);
	param->get_boolean("varyrec",varyrec,false,false);

	if (bstart > 0 && bend < gen){
		bkpt = true;
	}

	string s;
	param->get_string ("initfreq", s, true);
	io::println ("initfreq=  " + s ,2);
	vector<string> toks; 
	functions::tokenize (s,  toks, ",");
	for (int i = 0 ; i < toks.size(); i++)
		initfreq.push_back (atof (toks[i].c_str()));

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, seed);

	
	param->print_parameters();

}

void twopointsim::runmany ()  {
	for (double r1 = 0.01; r1  <= 0.5; r1 += 0.01) { 
		recrate=r1*1e-2;
		vector<double> f00(ensize);
		vector<double> f01(ensize);
		vector<double> f10(ensize);
		vector<double> f11(ensize);
		vector<int> d00(ensize);
		vector<int> d01(ensize);
		vector<int> d10(ensize);
		vector<int> d11(ensize);
		for (int i = 0 ; i < ensize; i++) {
			f00[i] = initfreq[0];
			f01[i] = initfreq[1];
			f10[i] = initfreq[2];
			f11[i] = initfreq[3];

		}
		double p[4];
		unsigned int draw[4];


		double d1 = 0 ;
		double d2 = 0;
		double rho1 = 0 ; 
		double rho2 = 0;
		int d2denom = 0;

	for (int t = 0 ; t < gen; t++){
		int n = popsize;
		if (bkpt) { 
			if ( t>= bstart && t<bend) {
				n = bsize;
			}
		}
		d1 = 0 ;
		d2 = 0;
		rho1 = 0 ; 
		rho2 = 0;
		d2denom = 0;
		for (int i = 0 ; i < ensize; i++) {
			double d = f00[i]*f11[i] - f01[i]*f10[i];
			double p1 = f11[i] + f10[i];
			double p2 = f11[i] + f01[i];
			double rho =  0;
			if ( p1 > 0 && p1 < 1 && p2 > 0  && p2 < 1) 
				rho = d / sqrt (p1*(1-p1)*p2*(1-p2));

			d1 += d;
			rho1 += rho;

			if (t==0) { 
			d00[i] = ((int)(popsize*initfreq[0]));
			d01[i] = ((int)(popsize*initfreq[1]));
			d10[i] = ((int)(popsize*initfreq[2]));
			d11[i] = ((int)(popsize*initfreq[3]));
			}

			bool ismono = (d01[i]+d11[i]==0 || d10[i]+d11[i]==0 || d00[i] + d01[i] == 0 || d00[i] + d10[i]==0);

			if (!ismono){
				d2 += d;
				rho2 += rho;
				d2denom++;
			}
			p[0] = f00[i] - recrate * d;
			p[1] = f01[i] + recrate * d;
			p[2] = f10[i] + recrate * d;
			p[3] = f11[i] - recrate * d;


			gsl_ran_multinomial (r, 4, n, p, draw);
			f00[i]=1.0*draw[0]/n;
			f01[i]=1.0*draw[1]/n;
			f10[i]=1.0*draw[2]/n;
			f11[i]=1.0*draw[3]/n;
			d00[i] =   draw[0];
			d01[i] = draw[1];
			d10[i] = draw[2];
			d11[i] = draw[3];
			
			if (io::debug>=1) {
				cout << draw[0] << "," << draw[1] << "," << draw[2] << "," << draw[3] << endl;
				cout << f00[i] << "," << f01[i] << "," << f10[i] << "," << f11[i] << endl;
			}
		}
		d1/=ensize;
		rho1/=ensize;
		if (d2denom>0) {
			d2/=d2denom;
			rho2 /=d2denom;
		}
		
	}
	cout <<  r1 <<  "\t" << d1 << "\t" << d2 << "\t" << rho1 << "\t" << rho2 << "\t" << d2denom << endl;
	}
}

void twopointsim::run ()  {
	vector<double> f00(ensize);
	vector<double> f01(ensize);
	vector<double> f10(ensize);
	vector<double> f11(ensize);
	vector<int> d00(ensize);
	vector<int> d01(ensize);
	vector<int> d10(ensize);
	vector<int> d11(ensize);
	for (int i = 0 ; i < ensize; i++) {
		f00[i] = initfreq[0];
		f01[i] = initfreq[1];
		f10[i] = initfreq[2];
		f11[i] = initfreq[3];

	}
	double p[4];
	unsigned int draw[4];

	for (int t = 0 ; t < gen; t++){
		int n = popsize;
		if (bkpt) { 
			if ( t>= bstart && t<bend) {
				n = bsize;
			}
		}
		double d1 = 0 ;
		double d2 = 0;
		double rho1 = 0 ; 
		double rho2 = 0;
		int d2denom = 0;
		for (int i = 0 ; i < ensize; i++) {
			double d = f00[i]*f11[i] - f01[i]*f10[i];
			double p1 = f11[i] + f10[i];
			double p2 = f11[i] + f01[i];
			double rho =  0;
			if ( p1 > 0 && p1 < 1 && p2 > 0  && p2 < 1) 
				rho = d / sqrt (p1*(1-p1)*p2*(1-p2));

			d1 += d;
			rho1 += rho;

			if (t==0) { 
			d00[i] = ((int)(popsize*initfreq[0]));
			d01[i] = ((int)(popsize*initfreq[1]));
			d10[i] = ((int)(popsize*initfreq[2]));
			d11[i] = ((int)(popsize*initfreq[3]));
			}

			bool ismono = (d01[i]+d11[i]==0 || d10[i]+d11[i]==0 || d00[i] + d01[i] == 0 || d00[i] + d10[i]==0);

			if (!ismono){
				d2 += d;
				rho2 += rho;
				d2denom++;
			}
			p[0] = f00[i] - recrate * d;
			p[1] = f01[i] + recrate * d;
			p[2] = f10[i] + recrate * d;
			p[3] = f11[i] - recrate * d;


			gsl_ran_multinomial (r, 4, n, p, draw);
			f00[i]=1.0*draw[0]/n;
			f01[i]=1.0*draw[1]/n;
			f10[i]=1.0*draw[2]/n;
			f11[i]=1.0*draw[3]/n;
			d00[i] =   draw[0];
			d01[i] = draw[1];
			d10[i] = draw[2];
			d11[i] = draw[3];
			
			if (io::debug>=1) {
				cout << draw[0] << "," << draw[1] << "," << draw[2] << "," << draw[3] << endl;
				cout << f00[i] << "," << f01[i] << "," << f10[i] << "," << f11[i] << endl;
			}
		}
		d1/=ensize;
		rho1/=ensize;
		if (d2denom>0) {
			d2/=d2denom;
			rho2 /=d2denom;
		}
		cout << "***\t" <<  t << "\t" << d1 << "\t" << d2 << "\t" << rho1 << "\t" << rho2 << "\t" << d2denom << endl;
		
	}
}

int main (int argc, char *argv[]) { 
	if (argc > 1) { 
		twopointsim ts = twopointsim (string(argv[1]));
		if (ts.varyrec) { 
			ts.runmany ();
		} else { 
			ts.run ();
		}
	}	
}

