#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "vectorfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "convertf.h"
#include "ram.h"
#include "ap.h"
#include "optimization.h"
#include "nag.h"
#include "nag_stdlib.h"
#include "nage04.h"
#include "nagx02.h"

using namespace alglib;

class kimurafit  { 
    public:
        void read_afs ( string filename , vector<double> &x, vector<double> &y)  ;
        double neglogl_model1 (double f, double t, vector<double> &par ) ;
        double neglogl_model2 (double f, double t, vector<double> &par ) ;
        double neglogl_model3 (double f, double t, vector<double> &par ) ;
        double neglogl_model4 (double f, double t, vector<double> &par ) ;

        double neglogl_model2 (double f, double t, double a, double b, double c, bool compute);

        double findmle_nm(vector<double> &par) ;
        double findmle_nm_internal(vector<double> &par) ;
        double findmle_grid ( vector<double> &par ) ; 

        double findpmle_t(double tau, vector<double> &par) ;
        double findpmle_f(double f, vector<double> &par) ;

        void read_input (string filename) ;
        void set_data ( vector<double> &tx, vector<double> &ty) ;
        void run_theory (double f, double tau) ;
        void simulate_theory () ;
        void simulate_theory_model1 () ;
        void simulate_theory_model2 () ;
        void simulate_theory_model3 () ;

        double findmle ( vector<double> &par)  ;
        void run ()  ;

        kimurafit (int col  = 1, double l = 0, double h = 1, double scale = 1); 
        kimurafit (int argc, char *argv[]) ;

        int col;
        int col1;
        int samplesize;
        double l;
        double h;
        int normalize;
        double scale ;
        int nmin ;
        int nmax;
        bool profile;
        vector<double> x;
        vector<double> y;
        vector<double> z;
        double intercept;
        vector<double> logposterior;
        vector<double> lpostf;
        vector<double> lpostt;

        vector<double> theoryy;
        vector<double> expectedy;
        vector<double> posty;
        vector<double> postmeany;
        data *d, *p;
        double f;
        double t;


        int seed;
		const gsl_rng_type * rng_T;
        gsl_rng * rng_r;

        int n ;
        double alpha;
        double tau;
        double a ;
        double sigma2;
        double post;
        double salpha;
        double stau;

        int proposal;
        bool accept;
        int mcmc_iters;
        int burnin;
        int printint;

        void mcmc ();
        void mcmc_init ();
        void mcmc_update ();
        double posterior ( double a, double alpha, double tau, double sigma2)  ;
        double lposterior ( double a, double alpha, double tau, double sigma2)  ;
        double lposterior ( double alpha, double tau)  ;
        double lposterior_model3 ( double a, double alpha, double tau) ;

        void simulate () ;

        string task;
        string init;
        bool fixnuisance;
        double postmeanalpha;
        double postmeantau;

        double tmin  ;
        double tmax ;
        double tinc ;
        double finc ;
        double fmin ;
        double fmax ;

        string optimizer;
        int optim_iters;
        int model;

        minlbfgsstate state;
        minlbfgsreport rep;
        //
        // Timer
        struct timeval now;
        long nowtime;

};


kimurafit::kimurafit (int col , double l, double h, double scale )  { 
    this->col = col;
    this->l = l ;
    this->h= h;
    this->scale = scale;
    this->nmin = -1;
    this->nmax = -1;
    tmin  = 0.01;
    tmax =  0.30;
    fmin = 0.01;
    fmax = 0.04;
    tinc = 0.01;
    finc = 0.01;

}

kimurafit::kimurafit (int argc, char *argv[]) {
	static const char *optString = "vhs";

	static const struct option longOpts[] = {
		{ "parameter", required_argument, NULL, 'p' },
		{ "verbose", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
        { "task", required_argument, NULL, 's'},
		{ NULL, no_argument, NULL, 0 }
	};

    d = new data (argc, argv, optString, longOpts);
	string pfile;
	d->get_string ("parameter",pfile,true);
    task = "mcmc";
    d->get_string ("task",  task, false);

	p = new data (pfile);
	p->print_parameters();
    p->get_int("col",col,1);
    p->get_int("col1",col1,-1);
    p->get_double("l",l,0.01);
    p->get_double("h",h,0.3);
    p->get_int ("samplesize",samplesize,-1);
    p->get_int ("debug", io::debug, 0);
    p->get_int ("normalize", normalize, 0);
    p->get_boolean ("profile", profile, false);
    
    string datafile ;
    p->get_string ("datafile",datafile, true);
    read_afs ( datafile, x, y);

    p->get_boolean ("fixnuisance", fixnuisance, false, false);
	p->get_int ("mcmc_iters", mcmc_iters, 100);
	p->get_int ("optim_iters", optim_iters, 100);
    optimizer = "nm";
    p->get_string ("optimizer", optimizer, false);
    p->get_int ("model", model, 1, false);

    p->get_double ("fmin",fmin,0.01,false);
    p->get_double ("fmax",fmax,0.04,false);
    p->get_double ("tmin",tmin,0.01,false);
    p->get_double ("tmax",tmax,0.30,false);
    p->get_double ("tinc",tinc,0.01,false);
    p->get_double ("finc",finc,0.005,false);

	p->get_int ("burnin", burnin, 0);
	p->get_int ("print_interval", printint, 1);
	if (burnin > mcmc_iters) {
		cerr << "Number of iterations must be > burnin" << endl;
		exit(1);
	}
	p->get_int ("seed",seed,1,false);
    p->get_int ("proposal", proposal, 1, false);
    init = "random";
    p->get_string ("init", init, false);

    p->get_double ("salpha", salpha, 0.1, false );
    p->get_double ("stau", stau, 0.1, false );

    gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	rng_r = gsl_rng_alloc (rng_T);
	gsl_rng_set (rng_r, seed);

    
}

void kimurafit:: run ()  {
    if (task.compare("simulate")==0)
        simulate ();
    else if (task.compare("mcmc")==0)
        mcmc();
    else if (task.compare ("theory")==0) 
        simulate_theory ();
    else if (task.compare ("mle") ==0) {
        vector<double> mle;
        findmle(mle);
    } else if (task.compare ("evaluate")==0){ 
        p->get_double ("a",a,a);
        p->get_double ("alpha",alpha,alpha);
        p->get_double ("tau",tau, tau);
        vector<double> par;
        double nll = neglogl_model3 (alpha, tau, par);

        run_theory (alpha,tau);
        cout << "f = " << alpha << "\tt = " << t << "\ta = " << par[0] << "\tnll = " << nll <<  endl;
        ofstream ofs;
        ofs.open ("eval.afs");
        for(int k = 0  ; k < theoryy.size(); k ++) { 
            double mu = par[0] * theoryy[k]; 
            double res = (y[k]-mu); res = res*res;
            res/= z[k]; res/=2;

            ofs << x[k] << "\t" <<  theoryy[k] << "\t" << res <<  endl;
        }
        ofs.close();

 
    } else {
        vector<double> mle;
        findmle(mle);
    }
}

void kimurafit::read_afs ( string filename , vector<double> &x, vector<double> &y)  {
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
    nmin = -1;
    nmax = -1;
    int n = 0 ;
    samplesize;
    
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		vector<string> toks;
		functions::tokenize(line.c_str(),toks," \t");
        x.push_back (atof(toks[0].c_str()));
        y.push_back (atof(toks[col].c_str()));
        if (col1 < toks.size() && col1 >= 0) {
            z.push_back (atof(toks[col1].c_str()));
        }
        if (x[n]>=l && nmin <0)
            nmin = n;
        if (x[n]>=h && nmax<0)
            nmax = n;
        n++;
	 
    }

    inp.close ();
    if (nmax  < 0)
        nmax = n - 1;
    samplesize = n  - 1 ;
    cout << "samplesize = " << samplesize << endl;

    double intercept = 0 ;
    for (int i  = nmin ;  i<=nmax;i++)
        intercept += y[i];
    intercept /= (nmax - nmin+1);


    if (normalize==1) { 
        double sum =  0 ;
        for (int i = 0;  i < y.size(); i++) { 
            if (x[i]>0 && x[i]<1)
                sum += y[i];
        }
       for (int i = 0;  i < y.size(); i++) { 
            if (sum>0 && x[i]>0 && x[i]<1)
                y[i]/=sum;
        } 
    }

    cout << "nmin = " << nmin << "\t nmax = " << nmax << endl;
    cout << "intercept =  " << intercept << endl;

}

void kimurafit::set_data ( vector<double> &tx, vector<double> &ty) { 

    x.resize(0);
    y.resize(0);
    for (int i  = 0 ; i < tx.size(); i++) {
        if (tx[i] >= l && nmin < 0)
            nmin = i;
        if (tx[i] >=h && nmax < 0 ) 
            nmax = i;     
        x.push_back (tx[i]);
        y.push_back (ty[i]);
    }
    if (nmax < 0 )
        nmax = tx.size()-1;
}

void kimurafit::run_theory (double f, double tau) { 
    ostringstream oss ;
    oss << "~np29/bin/kimtrans " << f <<" " << tau << ">tmp";

    ostringstream ss;
    ss << "tmp";

    system (oss.str().c_str());
    ifstream inp (ss.str());
    if (!inp.is_open()){
        cerr << "Error reading file theory"  <<endl;
        exit(1);
    }
    string line;
    theoryy.resize(0);
    double norm ;
    double total   = 0;
    int linenum = 0 ;
    vector<double> tmpx;
    vector<double> tmpy;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#') {
            if (linenum == 2) { 
                vector<string> toks;
                functions::tokenize(line.c_str(),toks," \t");
                double f0 = atof (toks[toks.size()-2].c_str());
                double f1 = atof (toks[toks.size()-1].c_str());
                if (f0<=0 && std::isnan(f1) ) f1=1;
                if (f1<=0 && std::isnan(f0) ) f0=1;
                norm = 1 -f0 -f1;
            }
            continue;
        }
		vector<string> toks;
		functions::tokenize(line.c_str(),toks," \t");
        tmpx.push_back (atof(toks[0].c_str()));
        double y= atof(toks[1].c_str());
        total += y;
        tmpy.push_back(y);
    }
    
    theoryy.resize(x.size());
//    samplesize =  0;
    if ( samplesize  > 0 ) {
        if (total > 0 )
            norm /= total;
        int denom = tmpy.size();
        if (io::debug >= 1)
            cout << "denom= " << denom << endl;
        for (int j =  0 ; j < tmpy.size(); j++){ 
            tmpy[j] /=denom;
        }
        for (int i =  0 ; i < x.size(); i++){ 
            double w = 0 ;
            int f = x[i]*samplesize;
            for (int j =  0 ; j < tmpx.size(); j++){ 
                double a = tmpx[j];
                w += pow(a,f) * pow(1-a,samplesize-f) * tmpy[j];
            }
            w *= mathfn::choose(samplesize,f);
            theoryy[i] = w;
        }
    } else  {
        total = 0 ;
        for (int i =  0, j = 0  ; i < x.size(); i++){ 
            while ( tmpx[j] < x[i]) j++;
            if (tmpx[j] == x[i])  {
                theoryy[i] = tmpy[j];
            } else {
                double w = (x[i]-tmpx[j-1])/(tmpx[j]-tmpx[j-1]);
                theoryy[i] = w*(tmpy[j]-tmpy[j-1]) + tmpy[j-1];
            }
            total += theoryy[i];
        }
        if (total > 0 )
             norm /=total;
        for (int i =  0 ; i < x.size(); i++)
            theoryy[i]  *= norm;

    }
    if (normalize==1) { 
        double sum = 0 ;
        for (int i = 0 ; i < theoryy.size(); i++)  {
            if (x[i] > 0 && x[i] < 1)
                sum += theoryy[i];
        }
        if (sum > 0 ) {
            for (int i = 0 ; i < theoryy.size(); i++)  {
                if (x[i] > 0 && x[i] < 1)
                    theoryy[i] /= sum;
            }
        }
    }
   
    if (io::debug>=1) {
        cout << "**Continuous AFS\t"<<f<<"\t"<<t<<"\n";
        for (int i = 0 ; i < tmpy.size(); i++)  {
            cout << tmpx[i] << "\t" << tmpy[i] << endl;
        }
        cout << "**Frequency\tObserved\tExpected\n";
        for (int i = 0 ; i < theoryy.size(); i++)  {
            cout << x[i] << "\t" << y[i] << "\t" <<  theoryy[i] << endl;
        }
    }


    inp.close ();
}

void kimurafit::simulate_theory () { 
    switch (model) { 
        case 1:
            simulate_theory_model1 (); break;
        case 2:
            simulate_theory_model2 (); break;
        case 3:
            simulate_theory_model3 (); break;
        default:
            break;
    }
}

void kimurafit::simulate_theory_model1 () { 
    double alpha  = 0.02;
    double tau = 0.1;
    p->get_double("alpha",alpha,0.02);
    p->get_double("tau",tau,0.3);

    cout << alpha << "\t" << tau << endl;
    run_theory (alpha, tau);

    ofstream ofs ("theory.afs");
    for(int k = 0  ; k < theoryy.size(); k ++) { 
        ofs << x[k] << "\t" <<  theoryy[k] << endl;
    }
    ofs.close ();

}

void kimurafit::simulate_theory_model2 () { 
    double alpha  = 0.02;
    double tau = 0.1;
    double tau1 = 0.3;
    p->get_double("alpha",alpha,0.02);
    p->get_double("tau",tau,0.3);
    p->get_double("tau1",tau1,0.3);

    cout << alpha << "\t" << tau << "\t" << tau1 <<  endl;
    run_theory (alpha, tau);

    ofstream ofs ("theory.afs");
    for(int k = 0  ; k < theoryy.size(); k ++) { 
        ofs << x[k] << "\t" <<  theoryy[k] << endl;
    }
    ofs.close ();

    ofs.open ("theory1.afs");
    for(int k = 0  ; k < theoryy.size(); k ++) { 
        double y =  (exp(-tau1)+exp(-3*tau1)) - exp(-3*tau1)*(samplesize*x[k]+1)/(samplesize+2);
        y/=(samplesize+1);
        y += theoryy[k];
        ofs << x[k] << "\t" <<  y << endl;
    }
    ofs.close ();

}

void kimurafit::simulate_theory_model3 () { 
    double alpha  = 0.02;
    double tau = 0.1;
    double theta ;
    p->get_double("alpha",alpha,0.02);
    p->get_double("tau",tau,0.1);
    p->get_double ("theta", theta,10, true);

    cout << alpha << "\t" << tau << endl;
    run_theory (alpha, tau);

    ofstream ofs ("theory.afs");
    for(int k = 0  ; k < theoryy.size(); k ++) { 
        double mean = theoryy[k] * theta; 
        double sd = z[k];  sd=pow(sd,0.5);
        int c = (gsl_ran_gaussian (rng_r, sd)  + mean);
        ofs << x[k] << "\t" <<  c << "\t" << z[k] << endl;
    }
    ofs.close ();

}

void kimurafit::simulate () { 
    double a = 5e5;
    double alpha  = 0.02;
    double tau = 0.1;
    double var = 5e4;

    run_theory (alpha, tau);
    for (int i = 0  ; i < 1000; i++) {
        double x = (1.*i)/1000;
        double y = a * theoryy[i] + gsl_ran_gaussian (rng_r, pow(var,0.5));
        cout << x << "\t" << y << endl;
    }
}

double kimurafit::lposterior ( double alpha, double tau)  {
   double sum =  0 ;
    if (alpha <0 || alpha > 1)
        return -1./0 ;
    if (tau < 0 )
        samplesize =  0;
    samplesize =  0;
        return -1./0 ;
        
    run_theory (alpha, tau);
    double logp  = 0 ;
    double y2 = 0  ; double k2 =  0; 
    double yk2 = 0;

    int n =  0;
    for (int i  = nmin ;  i<=nmax;i++) { 
        y2 += y[i]*y[i]; 
        k2 += theoryy[i]*theoryy[i];
        yk2 += y[i] * theoryy[i];
        n++;
    }
    if (n > 0 ) { 
        y2 /= n; 
        k2 /= n;
        yk2 /= n;
    }
    double c = y2 - (yk2*yk2)/k2;
    c = c*0.5*n;
    logp = gsl_sf_lngamma (0.5*(n+1)) - 0.5*(n+1)*log(c) -0.5 * log (k2);

    return logp;

}




double kimurafit::lposterior_model3 ( double a, double alpha, double tau)  {
   double sum =  0 ;
    if (alpha <0 || alpha > 1)
        return -1./0 ;
    if (tau < 0 )
        return -1./0 ;
        
    run_theory (alpha, tau);
    double logp  = 0 ;

    int n =  0;
    for (int i  = nmin ;  i<=nmax;i++) { 
        double mu = a*theoryy[i];
        double var  = z[i];
        double tmp = (y[i]-mu); tmp = tmp * tmp; tmp =-tmp;
        tmp/=var; tmp/=2;
        logp += tmp;
    }
    return (logp);
 
}

double kimurafit::lposterior ( double a, double alpha, double tau, double sigma2)  {
   double sum =  0 ;
    if (alpha <0 || alpha > 1)
        return -1./0 ;
    if (tau < 0 )
        return -1./0 ;
    return lposterior (alpha, tau);
        
    run_theory (alpha, tau);
    double logprior  =  -log(sigma2);
    double logl   = 0 ;
    for (int i  = nmin ;  i<=nmax;i++) { 
        logl = logl -0.5*pow((y[i]-a*theoryy[i]),2)/sigma2 - 0.5 * log(sigma2);
    }
    double logpost = logl + logprior;
    return logpost;
}

double kimurafit::posterior ( double a, double alpha, double tau, double sigma2)  {
    if (alpha <0 || alpha > 1)
        return 0 ;
    if (tau < 0 )
        return 0 ;
    double y =  lposterior(a,alpha,tau, sigma2);
    return exp(y);
}

/*
double neglogd (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(gsl_vector_get(v,0)); f = (gp->fmin+gp->fmax*f)/(1+f);
    double t=exp(gsl_vector_get(v,1));
         
    gp->run_theory (f,t);
    double m1 = 0 ; double var = 0 ; double logl =  0;
    double mleintercept = 0 ;
    gp->findmarginal ( mleintercept, logl, m1, var);
    logl = -logl;

    return (logl);
}

double neglogd1 (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(gsl_vector_get(v,0)); f = (gp->fmin+gp->fmax*f)/(1+f);
    double t=exp(gsl_vector_get(v,1));

    gp->run_theory (f,t);

    double theta = 0 ;
    double theta1 =  0 ;  double theta2 = 0;
    for (int i  = gp->nmin ;  i<=gp->nmax;i++) {
        double a = gp->y[i];
        double b = gp->theoryy[i];
        theta1 += a*b;
        theta2 += b*b;
    }
    theta = theta1/theta2;


    if (io::debug >=1)  {
        cout << "f = "<<gp->f <<"\tt="<<gp->t<<endl;
    }

    for( int i = gp->nmin ; i <= gp->nmax; i++){
        double mu = theta*gp->theoryy[i];
        double var  = gp->z[i];
        double tmp = (gp->y[i]-mu) ; tmp = tmp*tmp ;
        tmp /= var; tmp /= 2;
        val += tmp;
    }
    return (val);
}
*/

double gsl_neglogl_model1 (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(gsl_vector_get(v,0)); f = (gp->fmin+gp->fmax*f)/(1+f);

    double t=exp(gsl_vector_get(v,1));
    t = (gp->tmin+gp->tmax*t)/(1+t);

    vector<double> par;
    return gp->neglogl_model1 ( f, t, par);
}

double gsl_neglogl_model1 (double x, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(x); f = (gp->fmin+gp->fmax*f)/(1+f);
    double t=gp->tau;

    vector<double> par;
    return gp->neglogl_model1 ( f, t, par);
}

double gsl_neglogl_model2 (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;

    double f=exp(gsl_vector_get(v,0)); 
    f = (gp->fmin+gp->fmax*f)/(1+f);

    double t=exp(gsl_vector_get(v,1));
    t = (gp->tmin+gp->tmax*t)/(1+t);

    vector<double> par;
    return gp->neglogl_model2 ( f, t, par);
}

double gsl_neglogl_model3 (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(gsl_vector_get(v,0)); f = (gp->fmin+gp->fmax*f)/(1+f);

    double t=exp(gsl_vector_get(v,1));
    t = (gp->tmin+gp->tmax*t)/(1+t);

    vector<double> par;
    return gp->neglogl_model3 ( f, t, par);
}


double gsl_neglogl_model3 (double x, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(x); f = (gp->fmin+gp->fmax*f)/(1+f);

    double t=gp->tau;

    vector<double> par;
    return gp->neglogl_model3 ( f, t, par);
}

double gsl_neglogl_model3_f (double x, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double t=exp(x); t = (gp->tmin+gp->tmax*t)/(1+t);

    double f=gp->alpha;

    vector<double> par;
    return gp->neglogl_model3 ( f, t, par);
}

double gsl_neglogl_model4 (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(gsl_vector_get(v,0)); f = (gp->fmin+gp->fmax*f)/(1+f);

    double t=exp(gsl_vector_get(v,1));
    t = (gp->tmin+gp->tmax*t)/(1+t);

    vector<double> par;
    return gp->neglogl_model4 ( f, t, par);
}


double gsl_neglogl_model4 (double x, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(x); f = (gp->fmin+gp->fmax*f)/(1+f);
    double t=gp->tau;

    vector<double> par;
    return gp->neglogl_model4 ( f, t, par);
}

double kimurafit::neglogl_model1 (double f, double t, vector<double> &par ) {
    run_theory (f,t);
    double a = 0; double  b = 0;
    int np = 0 ;
    double c2 = 0; double g2 = 0 ;
    double l = 0 ;
    for (int i = nmin ; i<=nmax; i++) { 
        a += y[i]*theoryy[i]; b+= theoryy[i]*theoryy[i];
        c2 += y[i]*y[i]; g2 += theoryy[i]*theoryy[i];
        np++;
    }
    c2/=np; g2/=np;

    double mleintercept = 0;

    if (b>0)
        mleintercept = a/b;

    for (int i =nmin ; i<=nmax; i++) {
        double tmp = y[i] - mleintercept * theoryy[i];
        l += tmp*tmp;
    }
    l /= np;
    double var = l ;
    double nl = 0.5*np*(1+log(l));


    double sigma2 = 1;
    double m1 =  -0.5*np*(c2-pow(mleintercept,2)*g2)/sigma2;
    m1 -= 0.5*(log(np)+log(g2));
    m1 -= 0.5*(np-1)*(log(sigma2));

    double marginal = c2 - pow(mleintercept,2) * g2;
    marginal  = -0.5*(np-1)*(log ( np) + log(marginal));
    marginal -= 0.5 * ( log(np)  + log(g2));
    marginal +=  gsl_sf_lngamma ((0.5*(np-1)));

    par.resize(4);
    par[0] = mleintercept;
    par[1] = m1;
    par[2] = var;
    par[3] = marginal;

    return nl;
}

double kimurafit::neglogl_model2 (double f, double t, vector<double> &par ) {
    run_theory (f,t);

    double check  = 0;
    for( int i = nmin ; i <= nmax; i++){
        check += theoryy[i];
    }
    int ncol  = 3;
    if (check == 0 )  {
        ncol = 2;
    } 

    int nrow = nmax - nmin + 1; 
    gsl_matrix *a = gsl_matrix_alloc ( ncol, ncol);
    gsl_vector *b = gsl_vector_alloc ( ncol);
    for (int  j = 0 ; j < ncol; j++){
        for (int k = 0 ; k < ncol; k++){
            gsl_matrix_set (a, j, k , 0 );
        }
        gsl_vector_set (b, j , 0);
    }

    double tmp[3];
    for( int i = nmin ; i <= nmax; i++){
        double var = 1/sqrt(z[i]);

        if (check == 0) { 
             tmp[0] = var; tmp[1] = -x[i]*var;
            for (int  j = 0 ; j < 2; j++){
                for (int k = 0 ; k < 2; k++){
                    double prev = gsl_matrix_get (a,j,k); prev += tmp[j]*tmp[k];
                    gsl_matrix_set ( a,j,k, prev);
                }
                double prev = gsl_vector_get (b,j);
                prev += y[i]*var*tmp[j];
                gsl_vector_set (b,j,prev);
            }

        } else { 
             tmp[0] = var; tmp[1] =theoryy[i]*var; tmp[2] = -x[i]*var;
            for (int  j = 0 ; j < 3; j++){
                for (int k = 0 ; k < 3; k++){
                    double prev = gsl_matrix_get (a,j,k); prev += tmp[j]*tmp[k];
                    gsl_matrix_set ( a,j,k, prev);
                }
                double prev = gsl_vector_get (b,j);
                prev += y[i]*var*tmp[j];
                gsl_vector_set (b,j,prev);
            }
        }
    }

    if (io::debug >= 2) { 
        cout << "check = " << check << endl;
        printf ("a=\n");
        gsl_matrix_fprintf (stdout, a, "%g");
        printf ("b=\n");
        gsl_vector_fprintf (stdout, b, "%g");
    }
    int s;
    gsl_vector *x0=  gsl_vector_alloc (ncol);
    gsl_permutation * p = gsl_permutation_alloc (ncol);
    gsl_linalg_LU_decomp (a, p, &s);
    gsl_linalg_LU_solve (a, p, b, x0);


    if (io::debug >= 2) { 
        printf ("x = \n");
        gsl_vector_fprintf (stdout, x0, "%g");
    }

    // Change of ordering !!!
    par.resize(3);
    if (check == 0  ) {
        par[1] = gsl_vector_get (x0,0);
        par[0] = 0;
        par[2] = gsl_vector_get (x0,1);
    } else {
        par[1] = gsl_vector_get (x0,0);
        par[0] = gsl_vector_get (x0,1);
        par[2] = gsl_vector_get (x0,2);
    }
    gsl_permutation_free (p);
    gsl_vector_free (x0);
    gsl_vector_free (b);
    gsl_matrix_free (a);

    double val = 0 ;
    for( int i = nmin ; i <= nmax; i++){
        double mu = par[0]*theoryy[i];
        mu += (par[1]-par[2]*x[i]);
        double var  = z[i];
        double tmp = (y[i]-mu) ; tmp = tmp*tmp ;
        tmp /= var; tmp /= 2;
        val += tmp;
    }
    return val;
}

double kimurafit::neglogl_model3 (double f, double t, vector<double> &par ) {
    double val  = 0;
    int np = 0 ;

    run_theory (f,t);

    double theta = 0 ;
    double theta1 =  0 ;  double theta2 = 0;
    for (int i  = nmin ;  i<=nmax;i++) {
        double a = y[i];
        double b = theoryy[i];
        theta1 += a*b;
        theta2 += b*b;
    }
    theta = theta2>0 ?theta1/theta2:0;

    par.resize(1);
    par[0] = theta;

    if (io::debug >=1)  {
        cout << "f = "<<f <<"\tt="<<t<<endl;
    }

    for( int i = nmin ; i <= nmax; i++){
        double mu = theta*theoryy[i];
        double var  = z[i];
        double tmp = (y[i]-mu) ; tmp = tmp*tmp ;
        tmp /= var; tmp /= 2;
        val += tmp;
    }
    return (val);
}

double kimurafit::neglogl_model4 (double f, double t, vector<double> &par ) {
    run_theory (f,t);
    int np = 0 ;
    double l = 0 ;
    for (int i =nmin ; i<=nmax; i++) {
        double tmp = y[i] - theoryy[i];
        l += tmp*tmp;
        np++;
    }
    l /= np;
    double var = l ;
    double nl = 0.5*np*(1+log(l));


    par.resize(1);
    return nl;
}


double kimurafit::neglogl_model2 (double f, double t, double a, double b, double c, bool compute) { 
    double val  = 0;
    int np = 0 ;
    if (compute)
        run_theory (f,t);

    double theta = c ;
    double theta1 =  0 ;  double theta2 = 0;
    for (int i  = nmin ;  i<=nmax;i++) {
        theta1 += y[i];
        theta2 += theoryy[i];
    }



    for( int i = nmin ; i <= nmax; i++){
        double mu = theta*theoryy[i];
        mu += (a-b*x[i]);
        double var  = z[i];
        double tmp = (y[i]-mu) ; tmp = tmp*tmp ;
        tmp /= var; tmp /= 2;
        val += tmp;
    }

    if (io::debug >=0)  {
        cout << "f = "<<f <<"\tt = "<<t<<"\t a ="<<a<<"\tb = "<< b << "\tc = " << c << "\tnll = " << val << endl;
    }
    return (val);
}

/*
double neglogd2 (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;
    double f=exp(gsl_vector_get(v,0)); f = (gp->fmin+gp->fmax*f)/(1+f);
    double t=exp(gsl_vector_get(v,1));
    double a=exp(gsl_vector_get(v,2));
    double b=exp(gsl_vector_get(v,3));
    double c=exp(gsl_vector_get(v,4));

    gp->run_theory (f,t);

    double theta = c ;
    double theta1 =  0 ;  double theta2 = 0;
    for (int i  = gp->nmin ;  i<=gp->nmax;i++) {
        theta1 += gp->y[i];
        theta2 += gp->theoryy[i];
    }


    if (io::debug >=1)  {
        cout << "f = "<<f <<"\tt="<<t<<endl;
    }

    for( int i = gp->nmin ; i <= gp->nmax; i++){
        double mu = theta*gp->theoryy[i];
        mu += (a-b*gp->x[i]);
        double var  = gp->z[i];
        double tmp = (gp->y[i]-mu) ; tmp = tmp*tmp ;
        tmp /= var; tmp /= 2;
        val += tmp;
    }
    return (val);
}

double neglogd3 (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    kimurafit *gp = (kimurafit *)param;

    double f=exp(gsl_vector_get(v,0)); 
    f = (gp->fmin+gp->fmax*f)/(1+f);

    double t=exp(gsl_vector_get(v,1));
    t = (gp->tmin+gp->tmax*t)/(1+t);

    vector<double> par;
    return gp->neglogl3 ( f, t, par);
}*/

double kimurafit::findmle ( vector<double> &par )  {
    if (optimizer.compare("grid")==0) { 
            findmle_grid(par);
    } else if (optimizer.compare("nm")==0){
            findmle_nm(par);
    }

    /*
    }else if (optimizer.compare ("lbfgs")==0) {
        findmle2_lbfgs(par);
    } else if (optimizer.compare ("nag") ==0 ) { 
        findmle2_nag(par);
    }*/
}

double kimurafit::findmle_grid ( vector<double> &par )  {
    double logl =  0;
    io::println ("Finding MLE",0);
        
    int iter =  0;
    bool flag = false;

    vector<double> fgrid;
    vector<double> tgrid;

    for (f = fmin ; f<=fmax; f+=finc) fgrid.push_back(f);
    for (t = tmin; t<=tmax; t+=tinc ) tgrid.push_back(t);

    double maxnll = 0;
    par.resize(4);
    double fmle = 0 ;
    double tmle = 0;
    int i = 0 ; int j = 0 ;
    vector<double> tprofile;
    vector<double> fprofile;

    for (f = fmin ; f<=fmax; f+=finc, i++, j = 0) { 
        for (t =tmin, j = 0 ; t<=tmax; t+=tinc , j++){ 
            vector<double> tmppar;
            double tmpnll = 0;
            if (model==1)
               tmpnll = neglogl_model1 (f, t, tmppar);
            else if (model == 2)
               tmpnll = neglogl_model2 (f, t, tmppar);
            else if (model == 3)
               tmpnll = neglogl_model3 (f, t, tmppar);
            else if (model == 4)
               tmpnll = neglogl_model4 (f, t, tmppar);


            if (model==1)
                cout << "f = " << f << "\tt=" << t << "\t" << tmppar[0] << "\t" << tmppar[1] << "\t" << tmppar[2] << "\t" << tmppar[3] << "\tnll =  "<<tmpnll <<  endl;
            else if (model==2)
                cout << "f = " << f << "\tt=" << t << "\t" << tmppar[0] << "\t" << tmppar[1] << "\t" << tmppar[2] << "\tnll =  "<<tmpnll <<  endl;
            else if (model==3)
                cout << "f = " << f << "\tt=" << t << "\t" << tmppar[0] << "\tnll =  "<<tmpnll <<  endl;
            else if (model==4)
                cout << "f = " << f << "\tt=" << t << "\tnll =  "<<tmpnll <<  endl;

            if (f <= fmin && profile)  {
                double tmppll = findpmle_t (t, tmppar);
                if (model==1)
                    cout << "profilell: t = " << t <<"\tf = " << tmppar[1]  << "\ta = " << tmppar[0] << "\tnll = " << tmppll << endl;
                else if (model==3)
                    cout << "profilell: t = " << t <<"\tf = " << tmppar[1]  << "\ta = " << tmppar[0] << "\tnll = " << tmppll << endl;
                else if (model==4)
                    cout << "profilell: t = " << t <<"\tf = " << tmppar[0] << "\tnll = " << tmppll << endl;
            }

            if (!std::isnan (tmpnll) && (flag == false || tmpnll  < maxnll)) {
                maxnll = tmpnll;
                if (model==1) { 
                    par[0] = tmppar[0]; 
                    par[1] = tmppar[1];
                    par[2] = tmppar[2];
                    par[3] = tmppar[3];
                } else if (model==2) { 
                    par[0] = tmppar[0]; 
                    par[1] = tmppar[1];
                    par[2] = tmppar[2];
                } else if (model==3) { 
                    par[0] = tmppar[0]; 
                }
                fmle = f;
                tmle = t;
                flag = true;
            }
            iter ++;
        }
    }
    double f0 = fmin ; double f1 = fmax;
    if (model==1)
        cout << "gridsearch MLE = " << fmle << "\t" << tmle << "\t" << par[0] << "\t" << par[1] << "\t" << par[2] <<  "\t" << par[3] << "\t" << maxnll <<  endl;
    else if (model==2)
        cout << "gridsearch MLE = " << fmle << "\t" << tmle << "\t" << par[0] << "\t" << par[1] << "\t" << par[2] <<  "\t" << maxnll <<  endl;
    else if (model==3)
        cout << "gridsearch MLE = " << fmle << "\t" << tmle << "\t" << par[0] <<  "\t" << maxnll <<  endl;

    fmin = fmle - finc; fmin = (fmin<0)?finc:fmin;  fmin = (fmin<=f0)?f0:fmin;
    fmax = fmle + finc; fmax = (fmax>=f1)?f1:fmax;
    tmin = tmle - tinc ; tmin = (tmin<0)?tinc:tmin; tmax = tmle + tinc;

    if (fmin==fmax) { 
        double nll = findpmle_f (fmin, par);
        if (model==1)
            cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << par[3] << "\t" << nll << endl;
        else if (model==2)
            cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << par[3] << "\t" << par[4] << "\t" << nll << endl;
        else if (model==3)
            cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << nll << endl;
    } else {
        double nll = findmle_nm_internal (par);
        if (model==1)
            cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << par[3] << "\t" << nll << endl;
        else if (model==2)
            cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << par[3] << "\t" << par[4] << "\t" << nll << endl;
        else if (model==3)
            cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << nll << endl;
    }
 

    run_theory (par[1],par[2]);
    ofstream ofs;
    ofs.open ("mle.afs");
    for(int k = 0  ; k < expectedy.size(); k ++) { 
        ofs << x[k] << "\t" <<  theoryy[k] << endl;
    }
    ofs.close();

    if( model==2) {
        ofs.open ("mle1.afs");
        for(int k = 0  ; k < expectedy.size(); k ++) { 
            ofs << this->x[k] << "\t" <<  (par[0]*theoryy[k] + par[3] - par[4]*this->x[k] )<< endl;
        }
        ofs.close ();
    }

}

/*
double kimurafit::findmarginal (double &mleintercept, double &mle, double &m1 , double &var) {
    double a = 0; double  b = 0;
    int np = 0 ;
    double c2 = 0; double g2 = 0 ;
    double l = 0 ;
    for (int i = nmin ; i<=nmax; i++) { 
        a += y[i]*theoryy[i]; b+= theoryy[i]*theoryy[i];
        c2 += y[i]*y[i]; g2 += theoryy[i]*theoryy[i];
        np++;
    }
    c2/=np; g2/=np;

    if (b>0)
        mleintercept = a/b;

    for (int i =nmin ; i<=nmax; i++) {
        double tmp = y[i] - mleintercept * theoryy[i];
        l += tmp*tmp;
    }
    l /= np;
    var = l ;
    mle = -0.5*np*(1+log(l));

    double sigma2 = 1;
    m1 =  -0.5*np*(c2-pow(mleintercept,2)*g2)/sigma2;
    m1 -= 0.5*(log(np)+log(g2));
    m1 -= 0.5*(np-1)*(log(sigma2));

    double marginal = c2 - pow(mleintercept,2) * g2;
    marginal  = -0.5*(np-1)*(log ( np) + log(marginal));
    marginal -= 0.5 * ( log(np)  + log(g2));
    marginal +=  gsl_sf_lngamma ((0.5*(np-1)));

    return marginal;
}

double kimurafit::findmle1(double &mleintercept, double &mlesigma2) {
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *steps;
    gsl_multimin_function my_func;

    my_func.n = 1;
    my_func.f = neglogd1;
    my_func.params = this;
    x = gsl_vector_alloc (1);
    steps = gsl_vector_alloc (1);
    gsl_vector_set (x, 0, intercept);
    gsl_vector_set (steps, 0, 100);
    T = gsl_multimin_fminimizer_nmsimplex;
    s = gsl_multimin_fminimizer_alloc (T, 1);
         
    gsl_multimin_fminimizer_set (s, &my_func, x, steps);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s),1e-3);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        if ( io::debug >= 1){
            printf ("Inner iter = %5d %.5f  %10.5f\n", iter,
                    gsl_vector_get (s->x, 0), 
                    s->fval);
         }

    }
    while (status == GSL_CONTINUE && iter < optim_iters);


    mleintercept = gsl_vector_get(s->x,0);
    mlesigma2 =  s->fval;

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);

    return (s->fval);
}

double kimurafit::findmle_nm(vector<double> &par) {
    par.resize(4);
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *steps;
    gsl_multimin_function my_func;

    my_func.n = 2;
    my_func.f = neglogd1;
    my_func.params = this;
    x = gsl_vector_alloc (2);
    steps = gsl_vector_alloc (2);

    double f = (fmin+fmax)*.5; double t = 0.1;
        
    io::println ("Finding MLE",0);

    gsl_vector_set (x, 0, log((f-fmin)/(fmax-f)));
    gsl_vector_set (x, 1, log(t));
    gsl_vector_set (steps, 0, 0.1);
    gsl_vector_set (steps, 1, 0.1);
    T = gsl_multimin_fminimizer_nmsimplex;
    s = gsl_multimin_fminimizer_alloc (T, 2);
         
    gsl_multimin_fminimizer_set (s, &my_func, x, steps);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s),1e-3);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        if ( io::debug >= 0){
            printf ("Inner iter = %5d %.5f %.5f %10.5f \n", iter,
                    gsl_vector_get (s->x, 0), 
                    gsl_vector_get (s->x, 1), 
                    s->fval);
         }

    }
    while (status == GSL_CONTINUE && iter < optim_iters);


    double mlef=exp(gsl_vector_get(s->x,0)); mlef = (fmin+fmax*mlef)/(1+mlef);
    double mlet  = gsl_vector_get(s->x,1);
    mlet = exp(mlet);

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);

    run_theory (mlef,mlet);
    expectedy.resize(theoryy.size());
    ofstream ofs ("mle.afs");
    for(int k = 0  ; k < expectedy.size(); k ++) { 
        ofs << this->x[k] << "\t" <<  theoryy[k] << endl;
    }
    ofs.close ();


    double mleintercept = 0 ;
    double logl =  0 ;
    double m1 = 0 ;double var = 0;
    double logmarginal  = findmarginal ( mleintercept, logl, m1, var);
    cout << "gridsearch:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t" << mleintercept << "\t" << logl << "\t" << m1 <<"\t"<< logmarginal << "\t" << var << endl;
    cout << "MLE = " << mlef << "\t" << mlet << "\t" << mleintercept << "\t" << var << endl;
    par[0] = mleintercept;
    par[1] = mlef;
    par[2] = mlet;
    par[3] = var;
}

void NAG_CALL _naglib_monit(const double fmin, const double fmax,
        const double sim[], const Integer n,
        const Integer ncall, const double serror,
        const double vratio, Nag_Comm *comm)
{
}

static void NAG_CALL _naglib_evaluate(Integer n,  const double x[], double *objf,
                                     Nag_Comm *comm)
{
    io::println ("In _naglib_evaluate" , 0);
    double fmin = reinterpret_cast<kimurafit*>(comm->p)->fmin;
    double fmax = reinterpret_cast<kimurafit*>(comm->p)->fmax;
    double f = exp(x[0]); f = (fmin + f*fmax)/(1+f);
    double t = exp(x[1]);
    double a = exp(x[2]);
    double b = exp(x[3]);
    double c = exp(x[4]);

    *objf = reinterpret_cast<kimurafit*>(comm->p)->neglogl2(f,t,a,b,c);
}

double kimurafit::findmle2_nag (vector<double> &par) {
    return findmle2_nag_internal (par);
}

double kimurafit::findmle2_nag_internal (vector<double> &par) {
    par.resize (5);

    Integer       exit_status = 0;
    Nag_Boolean   print;
    Integer       monitoring, maxcal = 200, n = 5;
    Nag_BoundType bound;
    Nag_E04_Opt   options;
    double        *bl = 0, *bu = 0, *g = 0, objf, *x = 0;
    NagError      fail;
    double      f, tolf, tolx;
    Nag_Comm *ng = NAG_ALLOC(1, Nag_Comm);
    ng->p = this;

    INIT_FAIL(fail);
    n = 5;
    if (n >= 1)
    {
        if (!(x = NAG_ALLOC(n, double)) ||
                !(g = NAG_ALLOC(n, double)) ||
                !(bl = NAG_ALLOC(n, double)) ||
                !(bu = NAG_ALLOC(n, double)))
        {
            printf("Allocation failure\n");
            exit_status = -1;
            goto END;
        }
    }
    else
    {
        printf("Invalid n.\n");
        exit_status = 1;
        return exit_status;
    }
    tolf = sqrt(nag_machine_precision);
    tolx = sqrt(tolf);
    x[0] = 0.5*(fmin + fmax); x[0] = log((x[0]-fmin)/(fmax-x[0]));
    x[1] = 0.3; x[1] =log(x[1]);
    x[2] = 10;
    x[3] = 10;
    x[4] = 10;

    nag_opt_init(&options);
    strcpy(options.outfile, "stdout");
    //nag_machine_precision (x02ajc).
    options.optim_tol = 100*sqrt(nag_machine_precision);

    // Read remaining option values from file
    print = Nag_TRUE;
    // Set bounds on variables
    bound = Nag_Bounds;
    bl[0] = fmin;
    bu[0] = fmax;
    bl[1] = 0;
    bu[1] = 1;
    bl[2] = -1.0e10;
    bu[2] = 1.0e10;
    bl[3] = -1.0e10;
    bu[3] = 1.0e10;
    bl[4] = -1.0e10;
    bu[4] = 1.0e10;

    // Call optimization routine 
    // nag_opt_bounds_no_deriv (e04jbc), see above. 
    //nag_opt_bounds_no_deriv(n, _naglib_evaluate, bound, bl, bu, x, &objf,
    //        g, &options, ng, &fail);
    nag_opt_simplex_easy (n, x, &objf, tolf, tolx, _naglib_evaluate, _naglib_monit, maxcal, ng, &fail);
    if (fail.code != NE_NOERROR)
    {
        printf("The final function value is %12.4f\n", f);
        printf("at the point");
        for (int i = 1; i <= n; ++i)
        {
            printf(" %12.4f", x[i-1]);
        }
        printf("\n");
        printf("Error/Warning from nag_opt_bounds_no_deriv (e04jbc).\n%s\n",
                fail.message);
        if (fail.code != NW_COND_MIN)
            exit_status = 1;
    }

    nag_opt_free(&options, "all", &fail);
    if (fail.code != NE_NOERROR)
    {
        printf("Error from nag_opt_free (e04xzc).\n%s\n", fail.message);
        exit_status = 2;
    }

END:
    if (x) NAG_FREE(x);
    if (g) NAG_FREE(g);
    if (bl) NAG_FREE(bl);
    if (bu) NAG_FREE(bu);

    return 0;
} */

double kimurafit::findmle_nm(vector<double> &par) {
    double nll =  findmle_nm_internal (par);
    if (model==1)
        cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << par[3] << "\t" << nll << endl;
    else if (model==2)
        cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << par[3] << "\t" << par[4] << "\t" << nll << endl;
    else if (model==3)
        cout << "MLE  = "  << par[1] << "\t" << par[2]  << "\t" << par[0] << "\t" << nll << endl;

    run_theory (par[1],par[2]);

    ofstream ofs;
    ofs.open ("mle.afs");
    for(int k = 0  ; k < expectedy.size(); k ++) { 
        ofs << x[k] << "\t" <<  theoryy[k] << endl;
    }
    ofs.close();

    if (model==2) {
        ofs.open ("mle1.afs");
        for(int k = 0  ; k < expectedy.size(); k ++) { 
            ofs << this->x[k] << "\t" <<  (par[0]*theoryy[k] + par[3] - par[4]*this->x[k] )<< endl;
        }
        ofs.close ();
    }
}


double kimurafit::findmle_nm_internal(vector<double> &par) {
    par.resize(2);
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *steps;
    gsl_multimin_function my_func;

    my_func.n = 2;
    if (model==1)
        my_func.f = gsl_neglogl_model1;
    else if (model==2)
        my_func.f = gsl_neglogl_model2;
    else if (model==3)
        my_func.f = gsl_neglogl_model3;

    my_func.params = this;
    x = gsl_vector_alloc (2);
    steps = gsl_vector_alloc (2);

    double f = (fmin+fmax)*.5; 
    double t = (tmin+tmax)*.5;
    cout << "fmin = " << fmin  << "\tfmax = " << fmax << endl;
    cout << "tmin = " << tmin  << "\ttmax = " << tmax << endl;

    gsl_vector_set (x, 0, log(f-fmin)-log(fmax-f));
    gsl_vector_set (x, 1, log(t-tmin)-log(tmax-t));

    gsl_vector_set (steps, 0, 10);
    gsl_vector_set (steps, 1, 0.1);

    T = gsl_multimin_fminimizer_nmsimplex;
    s = gsl_multimin_fminimizer_alloc (T, 2);
         
    gsl_multimin_fminimizer_set (s, &my_func, x, steps);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s),1e-3);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        if ( io::debug >= 0){
            printf ("Inner iter = %5d %.5f %.5f %10.5f \n", iter,
                    gsl_vector_get (s->x, 0), 
                    gsl_vector_get (s->x, 1), 
                    s->fval);
         }

    }
    while (status == GSL_CONTINUE && iter < optim_iters);


    double mlef=exp(gsl_vector_get(s->x,0)); mlef = (fmin+fmax*mlef)/(1+mlef);
    double mlet  = exp(gsl_vector_get(s->x,1)); mlet = (tmin + tmax*mlet)/(1+mlet);
    double nll = s->fval;


    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);

    run_theory (mlef,mlet);
    expectedy.resize(theoryy.size());
    
    if (model==1) {
        vector<double> tmppar;
        double logl = neglogl_model1 (mlef, mlet, tmppar);
        double mleintercept = tmppar[0];
        double m1 = tmppar[1];
        double var = tmppar[2];
        double logmarginal  = tmppar[3];

        cout << "nm:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t" << mleintercept << "\t" << logl << "\t" << m1 <<"\t"<< logmarginal << "\t" << var << endl;
        par.resize(4);
        par[0] = mleintercept;
        par[1] = mlef;
        par[2] = mlet;
        par[3] = var;

    } else if (model ==2 ){
        vector<double> tmppar;
        neglogl_model2 (mlef, mlet, tmppar);
        double mlea = tmppar[1]; double mleb = tmppar[2]; double mlec = tmppar[0];

        cout << "a = " << mlea << "\tb = " << mleb << endl;
        cout << "nm:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t" << mlec << "\t" << nll << "\t" << mlea <<"\t"<< mleb << "\t" << nll << endl;

        par.resize(5);
        par[0] = tmppar[0];
        par[1] = mlef;
        par[2] = mlet;
        par[3] = tmppar[1];
        par[4] = tmppar[2];
    } else if (model == 3)  {
        vector<double> tmppar;
        neglogl_model3 (mlef, mlet, tmppar);
        cout << "nm:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t"<< tmppar[0] << "\t" << nll << endl;

        par.resize(3);
        par[0] = tmppar[0];
        par[1] = mlef;
        par[2] = mlet;
    }

    return nll;
}

double kimurafit::findpmle_t(double tau, vector<double> &par) {
    size_t iter = 0;
    int status;

    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function my_func;
    this->tau = tau;
    if (model==1)
       my_func.function = gsl_neglogl_model1;
    else if (model==3)
       my_func.function = gsl_neglogl_model3;
    else if (model==4)
       my_func.function = gsl_neglogl_model4;

    my_func.params = this;
    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);

    double f = (fmin+fmax)*.5; 
    double a = -100 ; double b = 100; 
    double m =log(f-fmin)-log(fmax-f); 
    double m_expected = M_PI;
    cout << fmin << "\t" << fmax << endl;
    cout << gsl_neglogl_model3 (a, this) << "\t" << gsl_neglogl_model3(b,this) <<"\t" << gsl_neglogl_model3(m,this) << endl;

    gsl_min_fminimizer_set (s, &my_func, m, a, b);

    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);
        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);
        status 
                     = gsl_min_test_interval (a, b, 0.001, 0.0);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        if ( io::debug >= 1){
            printf ("Inner iter = %5d [%.5f %.5f] %10.5f %10.5f %10.5f \n", iter,
                            a,b,m,m-m_expected, b-a);
         }

    }
    while (status == GSL_CONTINUE && iter < optim_iters);


    double mlef=gsl_min_fminimizer_x_minimum (s);
    mlef=exp(mlef); mlef = (fmin+fmax*mlef)/(1+mlef);
    double nll=gsl_min_fminimizer_f_minimum(s);
    double mlet = tau;

    gsl_min_fminimizer_free (s);

    run_theory (mlef,t);
    expectedy.resize(theoryy.size());

    if (model == 1)  {
        vector<double> tmppar;
        neglogl_model1 (mlef, mlet, tmppar);
        cout << "profilell:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t"<< tmppar[0] << "\t" << nll << endl;

        par.resize(3);
        par[0] = tmppar[0];
        par[1] = mlef;
        par[2] = mlet;
    
    } else if (model == 3)  {
        vector<double> tmppar;
        neglogl_model3 (mlef, mlet, tmppar);
        cout << "profilell:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t"<< tmppar[0] << "\t" << nll << endl;

        par.resize(3);
        par[0] = tmppar[0];
        par[1] = mlef;
        par[2] = mlet;
    } else if (model == 4) { 
        vector<double> tmppar;
        neglogl_model4 (mlef, mlet, tmppar);
        cout << "profilell:iter = " << iter << "\t" <<mlef << "\t" << mlet << "\t" << nll << endl;
        par.resize(2);
        par[0] = mlef;
        par[1] = mlet;
    }

    return nll;
}


double kimurafit::findpmle_f (double f, vector<double> &par) {
    size_t iter = 0;
    int status;

    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function my_func;
    this->alpha = f;
    if (model==3)
       my_func.function = gsl_neglogl_model3_f;

    my_func.params = this;
    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);

    double t = (tmin+tmax)*.5; 
    double a = -100 ; double b = 100; 
    double m =log(t-tmin)-log(tmax-t); 
    double m_expected = M_PI;
    cout << gsl_neglogl_model3_f (a, this) << "\t" << gsl_neglogl_model3_f(b,this) <<"\t" << gsl_neglogl_model3_f(m,this) << endl;

    gsl_min_fminimizer_set (s, &my_func, m, a, b);

    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);
        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);
        status 
                     = gsl_min_test_interval (a, b, 0.001, 0.0);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        if ( io::debug >= 0){
            printf ("Inner iter = %5d [%.5f %.5f] %10.5f %10.5f %10.5f \n", iter,
                            a,b,m,m-m_expected, b-a);
         }

    }
    while (status == GSL_CONTINUE && iter < optim_iters);


    double mlet=gsl_min_fminimizer_x_minimum (s);
    mlet=exp(mlet); mlet = (tmin+tmax*mlet)/(1+mlet);
    double nll=gsl_min_fminimizer_f_minimum(s);
    double mlef = alpha;

    gsl_min_fminimizer_free (s);

    run_theory (mlef,mlet);
    expectedy.resize(theoryy.size());

    if (model == 1)  {
        vector<double> tmppar;
        neglogl_model1 (mlef, mlet, tmppar);
        cout << "profilell:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t"<< tmppar[0] << "\t" << nll << endl;

        par.resize(3);
        par[0] = tmppar[0];
        par[1] = mlef;
        par[2] = mlet;
    
    } else if (model == 3)  {
        vector<double> tmppar;
        neglogl_model3 (mlef, mlet, tmppar);
        cout << "profilell:iter = " << iter << "\t" <<mlef << "\t" << mlet <<"\t"<< tmppar[0] << "\t" << nll << endl;

        par.resize(3);
        par[0] = tmppar[0];
        par[1] = mlef;
        par[2] = mlet;
    } else if (model == 4) { 
        vector<double> tmppar;
        neglogl_model4 (mlef, mlet, tmppar);
        cout << "profilell:iter = " << iter << "\t" <<mlef << "\t" << mlet << "\t" << nll << endl;
        par.resize(2);
        par[0] = mlef;
        par[1] = mlet;
    }

    return nll;
}
/*
static void _alglib_evaluate (const real_1d_array &x, double &func, void *ptr)  {
    func = reinterpret_cast<kimurafit*> (ptr)->alglib_evaluate (x);
}

double kimurafit::alglib_evaluate (const real_1d_array &x ) { 

    double f=exp(x[0]); f = (fmin+fmax*f)/(1+f);
    double t=exp(x[1]);
    double a=exp(x[2]);
    double b=exp(x[3]);
    double c=exp(x[4]);

    return neglogl2 (f, t, a, b, c);
}

static void _alglib_progress (const real_1d_array &x, double func, void *ptr) { 
    reinterpret_cast<kimurafit*> (ptr)->alglib_progress (x, func );
}


void kimurafit::alglib_progress( const real_1d_array &x, double fx)  {
    alglib_impl::minlbfgsstate *p = state.c_ptr();
    int k = p->repiterationscount;
    gettimeofday(&now, NULL);
    long iterationtime = now.tv_sec * UMILLION + now.tv_usec;
    long elapsed = iterationtime - nowtime;
    elapsed/=1e6;
    nowtime = iterationtime;
    int n = 5;

    cout << "Iteration " << k << "\tfx = " << fx ;
    double xnorm = 0 ;
    for (int i  = 0 ; i < n ; i++) { 
        cout << "," << "x["<< i << "] = " << x[i] ;
        xnorm += x[i]*x[i];
    }
    cout << endl;
    xnorm = sqrt(xnorm);
    cout << "xnorm = " << xnorm << endl;
 
}

double kimurafit::findmle2_lbfgs(vector<double> &par) {
    real_1d_array x;
    x.setlength(5);
    double epsg = 0.0000000001;
    double epsf = 1e-4;
    double epsx = 0;
    ae_int_t maxits = 0;
    double f = (fmin+fmax)*.5; double t = 0.1;
    double a = 1e5; 
    double b = 1e5;
    double c = 1e5;
    x[0] = log((f-fmin)/(fmax-f));
    x[1] = log(t);
    x[2] = log(a);
    x[3] = log(b);
    x[4] = log(c);

    cout << "Initial = " ;
    cout << "x["<< 0 << "] = " << x[0] ;
    for (int i  = 1 ; i < 5 ; i++) { 
        cout << "," << "x["<< i << "] = " << x[i] ;
    }
    cout << endl;

    try {
        minlbfgscreatef(3, x, 1e-6, state);
        minlbfgssetcond(state, epsg, epsf, epsx, maxits);
        minlbfgssetxrep (state, 1);
        alglib::minlbfgsoptimize(state, _alglib_evaluate, _alglib_progress, this);
        minlbfgsresults(state, x, rep);
        double fx = state.f;
        io::println ("LBFGS terminated with " + tostring (rep.terminationtype), 0 );


    } catch(alglib::ap_error e)
    {
            printf("error msg: %s\n", e.msg.c_str());
    }
    
}
*/

void kimurafit::mcmc () { 
    io::println ("Starting MCMC",0);
    mcmc_init ();
    int denom  = 0 ;
    int nef = (mcmc_iters - burnin)/printint;
    postmeanalpha = 0 ;
    postmeantau = 0 ;
	for (int i  = 0 ; i  < mcmc_iters; i++) { 
        if ( i < burnin )  {
            cout << "mcmc:burnin = " << i << "\t" << a << "\t" << alpha << "\t" << tau << "\t" << sigma2 << "\t" << accept << "\t"<< post << endl;
        } else {
            cout << "mcmc:iter = " << i << "\t" << a << "\t" << alpha << "\t" << tau << "\t" << sigma2 << "\t" << accept<<"\t"<< post << endl;
            for (int j = 0 ; j < theoryy.size(); j ++) 
                posty[j] += theoryy[j];
            postmeanalpha += alpha;
            postmeantau += tau;
            denom ++;
        }
        mcmc_update ();
    }
    if (denom > 0 ) { 
        for (int i = 0 ; i < theoryy.size(); i++)
            posty[i] /= denom;
        postmeanalpha /= denom;
        postmeantau /= denom;

        run_theory (postmeanalpha, postmeantau);
        cout << "posterior mean = " << postmeanalpha << "\t" << postmeantau << endl;
        ofstream ofs1 ("post.afs");
        for(int k = 0  ; k < expectedy.size(); k ++) { 
            ofs1 << x[k] << "\t" <<  posty[k] << endl;
        }
        ofs1.close ();
        ofstream ofs2 ("postmean.afs");
        for(int k = 0  ; k < expectedy.size(); k ++) { 
            ofs2 << x[k] << "\t" <<  theoryy[k] << endl;
        }
        ofs2.close ();
    }
}

void kimurafit::mcmc_init ( ) {
    if (init.compare("random")==0){
        a = gsl_ran_flat (rng_r, 0 , 1);
        alpha = gsl_ran_flat (rng_r, 0 , 0.1);
        tau = gsl_ran_flat (rng_r, 0 , 0.05);
        sigma2 = 100;
    } else if (init.compare ("mle")==0){
        vector<double> mle;
        findmle (mle);
        a = mle[0]; 
        alpha = mle[1]; tau = mle[2];
        sigma2 = mle[3];
    } else if (init.compare ("file")==0){
        a = gsl_ran_flat (rng_r, 0 , 1);
        alpha = gsl_ran_flat (rng_r, 0 , 0.1);
        tau = gsl_ran_flat (rng_r, 0 , 0.05);
        sigma2 = 100;

        p->get_double ("a",a,a);
        p->get_double ("alpha",alpha,alpha);
        p->get_double ("tau",tau, tau);
        p->get_double ("sigma2", sigma2, sigma2);
    } else {
        a=205460; sigma2 = 5e4;
        alpha = 0.02188; tau = 0.0780846;
    }

    run_theory (alpha, tau);
    if (model==1)
        post =  lposterior (a, alpha, tau, sigma2);
    else
        post = lposterior_model3(a, alpha, tau);
    accept = true;
    posty.resize (theoryy.size());
    postmeany.resize(theoryy.size());
}



void kimurafit::mcmc_update () {
    int n =  0 ;
    double s = 0  ;
    s = 0; double mean =  0;
    run_theory (alpha, tau);
    for (int i  = nmin ;  i<=nmax;i++) { 
        n++;
        s +=  pow((theoryy[i]),2)/z[i];
        mean += y[i]*theoryy[i]/z[i];
    //    cout << y[i] << "\t" << a << "\t" << theoryy[i] << "\t" << s << endl;
    }
    if ( io::debug >= 1) { 
        cout << "s = " << s << endl;
    }
    mean /=s;
    if ( !fixnuisance){ 
        s = 1/s; s = pow(s,0.5);
        cout << "mean = " << mean << "\t s = " <<  s << endl; 
        a = gsl_ran_gaussian (rng_r, s) + mean;
    }

    double tmpalpha; 
    double tmptau;
    double r;
    switch (proposal) { 
        case 1:
            tmpalpha =  log(alpha/(1-alpha));
            tmpalpha  = gsl_ran_gaussian (rng_r, salpha) + tmpalpha;
            tmpalpha = 1/(1+exp(-tmpalpha));
            tmptau = log(tau);
            tmptau = gsl_ran_gaussian (rng_r, stau) + tmptau;
            tmptau = exp(tmptau);
        break;
        case 2:
            tmpalpha = alpha + gsl_ran_gaussian (rng_r, salpha);
            tmptau = tau + gsl_ran_gaussian (rng_r, stau);
        break;
    }
    
    double tmppost;
    if (model==1) {
        tmppost =  lposterior (a, tmpalpha, tmptau, sigma2);
    } else if (model==3) { 
        tmppost =  lposterior_model3 (a, tmpalpha, tmptau);
    }

    double r1;
    switch (proposal) { 
        case 1:
            if (std::isinf(tmppost) || std::isnan(tmppost))
                r = -1./0 ;
            else {
                r = log(mathfn::dsigmoid(alpha))+log(tmptau)-log(mathfn::dsigmoid(tmpalpha))-log(tau) +tmppost - post;
                r1= r;
                r = mathfn::min (r,0);
            }
        break;
        case 2:
            if (std::isinf(tmppost) ||  std::isnan(tmppost))
                r = -1./0 ;
            else {
                r = tmppost-post;
                r1= r;
                r = mathfn::min(r,0);
            }
        break;
    }

    if ( io::debug >= 0) 
        cout << "Proposed = "<<a << "\t" << tmpalpha << "\t" << tmptau << "\t" << sigma2 << "\t" << tmppost<<"\t"<<r1 <<"\t" << r << endl;
    double u  = gsl_ran_flat ( rng_r, 0 , 1);

    if (log(u) <= r)  {
        accept = true;
        alpha = tmpalpha; 
        tau = tmptau;  
        post = tmppost;
    } else  {
        accept = false;
    }
    

}

void kimurafit::read_input (string filename) { 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
    nmin = -1;
    nmax = -1;
    int n  = 0;
	while ( std::getline (inp, line)){
		linenum ++;
		char c = line[0];
		if (c=='#')
			continue;
		vector<string> toks;
		functions::tokenize(line.c_str(),toks," \t");
        x.push_back (atof(toks[0].c_str()));
        y.push_back (atof(toks[col].c_str()));
        if (x[n]>=l && nmin <0)
            nmin = n;
        if (x[n]>=h && nmax<0)
            nmax = n;
        n++;
	}
    inp.close ();
    if (io::debug >= 2){
        for (int i = 0 ; i < n; i++)
            cout << x[i] << "\t" << y[i] << endl;
    }
}


int main (int argc, char *argv[]) { 
    kimurafit k = kimurafit (argc, argv);
    k.run ();
}
