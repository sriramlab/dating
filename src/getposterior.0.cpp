#include "std.h"
#include "data.h"
#include "io.h"
#include "fileio.h"
#include "stringfn.h"
#include "vectorfn.h"
#include "mathfn.h"
#include "snpmap.h"
#include "convertf.h"
#include "arms/arms.h"

ofstream statsfs;

class getposterior {
    public:
        getposterior (int argc, char *argv[]);
        double findmle() ;
        void read_rolloff (string filename) ;
        void sample_ygf();
        void sample_tgf () ;
        void sample_alpha () ;
        void sample_sigma2 ();
        void sample_lambda () ;
        void sample_intercept () ;
        double mcmc_init () ;
        double mcmc();
        void stats (double *x, int n, double &mean, double &sd, double &l, double &u) ;

        double logd (double intercept, double lambda, double sigma2) ;
        double logd (double intercept, double lambda) ;
        double logd (double lambda) ;
        double logd1 (double t, double alpha) ;
        double logd2 (double y, double g, double alpha) ;
        int col;
        double l;
        double h;
        double s;
        vector<double> x;
        vector<double> y;
        int n;
        int nmin ;
        int nmax;
       	// Parameters passed to this program
		data *d;

		// Parameters specific to the data
		data *p;

        bool fitaffine;
        double affine; 

        double sigma2 ; 
        double intercept;
        double lambda;
        double alpha;
        double tgf;
        double ygf;
        double g[2];

        double mlelambda;
        double mleintercept;
        double mletgf ;
        double mlesigma2;

        double tail;

        vector<double> alphav;
        double *alphaprob;
        unsigned int *alphapick;


        int seed;
		const gsl_rng_type * rng_T;
        gsl_rng * rng_r;

        int mcmc_iters;
        int burnin;
        int printint;


        bool givenjack;
        double lambdahat;
        double lambdase;
};


getposterior::getposterior (int argc, char *argv[]) {
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
    p->get_int("col",col,2);
    p->get_double("l",l,0.02);
    p->get_double("h",h,1);
    string inpfile;
    p->get_string ("inp",inpfile,true);
    read_rolloff (inpfile);
    string alphafile;
    p->get_string ("alpha", alphafile, true);
    p->get_double ("glb", g[0], 25);
    p->get_double ("gub", g[1], 33);
    string statsfile;
    p->get_string ("output", statsfile, true);
    p->get_double ("tail", tail, 0.05);
    statsfs.open (statsfile.c_str());

    vector<string> *tmp = new vector<string>();
    fileio::read_vector (alphafile, tmp,0);

    alphav = vector<double> (tmp->size());
    alphaprob = new double [tmp->size()];
    alphapick = new unsigned int[tmp->size()];
    for (int i = 0 ; i  < tmp->size();i++){
        alphav[i] = atof((*tmp)[i].c_str());
        alphaprob[i] = 1.0/tmp->size();
        alphapick[i] = 0;
    }

	p->get_int ("debug", io::debug, 0);
	p->get_int ("mcmc_iters", mcmc_iters, 1200);
	p->get_int ("burnin", burnin, 200);
	p->get_int ("print_interval", printint, 1);
	if (burnin > mcmc_iters) {
		cerr << "Number of iterations must be > burnin" << endl;
		exit(1);
	}
	p->get_int ("seed",seed,1,false);


    p->get_boolean ("givenjack", givenjack, false, false);
    if (givenjack )  {
        p->get_double ("lambdahat", lambdahat, 1, true);
        p->get_double ("lambdase", lambdase, 1, true);
    }
    p->get_boolean ("fitaffine", fitaffine, false, false);

    gsl_rng_env_setup();
	rng_T = gsl_rng_default;
	rng_r = gsl_rng_alloc (rng_T);
	gsl_rng_set (rng_r, seed);

}



void getposterior::read_rolloff (string filename) { 
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

double getposterior::logd2 (double y, double g, double alpha) {
    double t = y/g;
    return logd1 (t,alpha);
}

double getposterior::logd1 (double t, double alpha) {
    double lambda = alpha * (log(t/alpha)+1);
    return logd (lambda);
}

double getposterior::logd (double lambda) {
    return logd (mleintercept, lambda, mlesigma2);
}

double getposterior::logd (double lambda, double intercept) {
    return logd (intercept, lambda, mlesigma2);
}

double getposterior::logd (double intercept, double lambda, double sigma2) {
    double val = 0;
    int np = 0 ;
    double affine = 0  ;
    if ( givenjack ) { 
        val = - 0.5 * pow  ( ( lambda - lambdahat ) , 2) / (lambdase*lambdase) - log  (lambdase);
        return (val);
    } else { 
    if (fitaffine)  {
        int nc = 0 ; 
        for( int i = nmin ; i <= nmax; i++){
            double tmp = y[i]-intercept*exp(-lambda*x[i]/100);
            affine += tmp;
            nc ++;
        }
        affine /= nc;
    } else 
        affine = 0 ;

    for( int i = nmin ; i <= nmax; i++){
        double tmp = y[i]-intercept*exp(-lambda*x[i]/100) - affine ;
        val += tmp*tmp;
        np++;
    }
//    val =  -0.5*np*log(val) ;
    val =  -0.5*val/sigma2 - 0.5 * np*log(sigma2) ;
    return (val);
    }
}

double neglogd (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    getposterior *gp = (getposterior *)param;
    double affine = 0  ;
    if (gp->fitaffine)  {
        int nc = 0 ; 
        for( int i = gp->nmin ; i <= gp->nmax; i++){
            double tmp = gp->y[i]-gsl_vector_get(v,0)*exp(-gsl_vector_get(v,1)*gp->x[i]/100);
            affine += tmp;
            nc ++;
        }
        affine /= nc;
    } else 
        affine = 0 ;

    for( int i = gp->nmin ; i <= gp->nmax; i++){
        double tmp = gp->y[i]-gsl_vector_get(v,0)*exp(-gsl_vector_get(v,1)*gp->x[i]/100) - affine ;
        val += tmp*tmp;
        np ++;
    }
    val = val *0.5/np;
    return (val);
}

double logllam (double lambda, void *param) { 
    getposterior *gp = (getposterior *)param;
    double val  = gp->logd (gp->intercept, lambda, gp->sigma2);
    return (val);

}

double loglint (double intercept, void *param){ 
    getposterior *gp = (getposterior *)param;
    double val  = gp->logd (intercept, gp->lambda, gp->sigma2);
    return (val);
}



double logltgf (double t, void *param) { 
    getposterior *gp = (getposterior *)param;

    double alpha = gp->alpha;
    vector<double> &av = gp->alphav;
    int an = av.size();
    vector<double> tmp (an);
    double ll  = 0;

    for (int i = 0 ; i  < an; i++) {
        double lambda = av[i] * log(t/av[i]+1);
        tmp[i] = gp->logd(gp->intercept,lambda, gp->sigma2);
    }
    double val  = vectorfn::lsumexp(tmp);
    return (val);

}


double getposterior::findmle() {
    size_t iter = 0;
    int status;


    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_vector *x;
    gsl_vector *steps;
    gsl_multimin_function my_func;

    my_func.n = 2;
    my_func.f = neglogd;
    my_func.params = this;
    x = gsl_vector_alloc (2);
    steps = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, 0.5);
    gsl_vector_set (x, 1, 100.0);
    gsl_vector_set (steps, 0, 0.01);
    gsl_vector_set (steps, 1, 10.0);
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

        if ( io::debug >= 1){
            printf ("%5d %.5f %.5f %10.5f\n", iter,
                    gsl_vector_get (s->x, 0), 
                    gsl_vector_get (s->x, 1), 
                    s->fval);
         }

    }
    while (status == GSL_CONTINUE && iter < 100);


    mleintercept = gsl_vector_get(s->x,0);
    mlelambda = gsl_vector_get(s->x,1);
    mlesigma2 =  s->fval;

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    return s->fval;
}


void getposterior::sample_sigma2 () {
    double val = 0;
    int np = 0 ;
    for( int i = nmin ; i <= nmax; i++){
        double tmp = y[i]-intercept*exp(-lambda*x[i]/100);
        val += tmp*tmp;
        np++;
    }
    val = 0.5*val;
    sigma2  = 1/gsl_ran_gamma (rng_r,0.5*np, 1/val);
}

void getposterior::sample_alpha () {
    double sum =  0;
    int  alphan= alphav.size();
    for (int i  = 0 ; i < alphan; i++) { 
        if (io::debug>=1)
            cout << tgf << "\t" << alphav[i] <<  "\t" << logd1(tgf,alphav[i]) << endl;
        alphaprob[i] = logd1 (tgf, alphav[i]);    
        sum += alphaprob[i];
        alphapick[i] = 0;
    }
    functions::normalize (alphaprob, alphan);
    for (int i  = 0 ; i < alphan; i++) { 
        alphaprob[i] = exp(alphaprob[i]);
    }

    if (io::debug>=1) { 
        for (int i = 0 ; i < alphav.size(); i++)
            cout << alphaprob[i] << "\t";
        cout << endl;
    }
    gsl_ran_multinomial (rng_r, alphav.size(), 1, alphaprob, alphapick);
    for (int i = 0 ; i  < alphav.size(); i++)  {
        if (alphapick[i]==1)
            alpha = alphav[i];
    }
}

void getposterior::sample_ygf() {
    ygf = gsl_ran_flat (rng_r, tgf*g[0], tgf*g[1]);
}

void getposterior::sample_tgf () {
	double xl = 0.1;
	double xr = 1e7;
	int ninit  = 30;
//	double xinit[7] = {1,10,100,1000,1e4,1e5,1e6};
	double xinit[30] =	{1,10,20,50,100,200,500,1e3,2e3,3e3,
				4e3,5e3,6e3,7e3,8e3,9e3,1e4,2e4,3e4,4e4,
				5e4,6e4,7e4,8e4,9e4,1e5,2.5e5,5e5,1e6,5e6};
	
	double xsamp[100];
	int npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
  	double xcent[10], qcent[10] = {5., 30., 70., 95.};
  	double convex = 1.;

	double xprev = tgf;
	double er = arms(xinit, ninit, &xl, &xr, logltgf, this,
			&convex,npoint, 1, &xprev, xsamp, nsamp, qcent, xcent, ncent, &neval);

	if (er > 0 ) { 
		cerr << "Error in sample_tgf: "  << er  << endl;
		exit (1);
	}

	if (io::debug >= 2){
		io::print ("ARMS:",1);
		for (int i =  0 ; i <  nsamp; i++)
			io::print (tostring(tostring(xsamp[i])) +",",1);
		io::print ("\n",1);
	}
	tgf = xsamp[0];
}

void getposterior::sample_lambda () {
    if (givenjack) { 
        lambda = gsl_ran_gaussian (rng_r, lambdase) + lambdahat ; 
    } else { 
	double xl = 0.1;
	double xr = 1e7;
	int ninit  = 30;
//	double xinit[7] = {1,10,100,1000,1e4,1e5,1e6};
	double xinit[30] =	{1,10,20,50,100,200,500,1e3,2e3,3e3,
				4e3,5e3,6e3,7e3,8e3,9e3,1e4,2e4,3e4,4e4,
				5e4,6e4,7e4,8e4,9e4,1e5,2.5e5,5e5,1e6,5e6};
	
	double xsamp[100];
	int npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
  	double xcent[10], qcent[10] = {5., 30., 70., 95.};
  	double convex = 1.;

	double xprev = lambda;
	double er = arms(xinit, ninit, &xl, &xr, logllam, this,
			&convex,npoint, 1, &xprev, xsamp, nsamp, qcent, xcent, ncent, &neval);

	if (er > 0 ) { 
		cerr << "Error in sample_lambda: "  << er  << endl;
		exit (1);
	}

	if (io::debug >= 2){
		io::print ("ARMS:",1);
		for (int i =  0 ; i <  nsamp; i++)
			io::print (tostring(tostring(xsamp[i])) +",",1);
		io::print ("\n",1);
	}
	lambda = xsamp[0];
    }
}

void getposterior::sample_intercept () {
	double xl = 0;
	double xr = 1;
	int ninit  = 30;
//	double xinit[7] = {1,10,100,1000,1e4,1e5,1e6};
	double xinit[30] =	{1,10,20,50,100,200,500,1e3,2e3,3e3,
				4e3,5e3,6e3,7e3,8e3,9e3,1e4,2e4,3e4,4e4,
				5e4,6e4,7e4,8e4,9e4,1e5,2.5e5,5e5,1e6,5e6};
    for (int i = 0  ;  i < 30; i ++){
        xinit[i] /= 5.1e6;
    }    
	double xsamp[100];
	int npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
  	double xcent[10], qcent[10] = {5., 30., 70., 95.};
  	double convex = 1.;

	double xprev = intercept;
	double er = arms(xinit, ninit, &xl, &xr, loglint, this,
			&convex,npoint, 1, &xprev, xsamp, nsamp, qcent, xcent, ncent, &neval);

	if (er > 0 ) { 
		cerr << "Error in sample_intercept: "  << er  << endl;
		exit (1);
	}

	if (io::debug >= 2){
		io::print ("ARMS:",1);
		for (int i =  0 ; i <  nsamp; i++)
			io::print (tostring(tostring(xsamp[i])) +",",1);
		io::print ("\n",1);
	}
	intercept = xsamp[0];
}

double getposterior::mcmc_init () {
    cout << "MLE = " << mleintercept << "\t" << mlelambda << endl;
    gsl_ran_multinomial (rng_r, alphav.size(), 1, alphaprob, alphapick);
    for (int i = 0 ; i < alphav.size(); i++)  {
        if (alphapick[i]==1)
            alpha = alphav[i];
    }

    tgf = alpha*(exp(mlelambda/alpha)-1);
    lambda = mlelambda;
    sigma2 = mlesigma2;
    io::println ("alphainit = " + tostring(alpha),0);
    io::println ("tgfinit = " + tostring(tgf),0);
    io::println ("lambdainit = " + tostring(lambda),0);
    io::println ("sigma2init  = " + tostring(sigma2),0);
}

double getposterior::mcmc ()  {
    io::println ("Starting MCMC",0);
    mcmc_init ();
    for (int i = 0  ; i < 100; i++){
        cout << loglint (i/100.0, this) << "\t";
    }
    cout <<endl;
    int ncollect = 3;
    double **collect = new double*[ncollect];
    int nef = (mcmc_iters - burnin)/printint;
    for (int  i  = 0 ; i  <ncollect;  i++)
        collect[i] = new double[nef];
    double* yl = new double[nef];
    double* yu = new double[nef];
    int j  = 0;
	for (int i  = 0 ; i  < mcmc_iters; i++) { 
        sample_intercept();
        sample_lambda();
        sample_sigma2 ();
        sample_tgf ();
        sample_ygf ();
        if ( i%printint == 0 )  {
			if ( i < burnin)
                cout << "burnin:" << i << "\t" << intercept << "\t" << lambda <<  "\t" << sigma2 << "\t" << tgf << "\t" << alpha << "\t" << ygf << endl;
            else {
                cout << "iter:" << i << "\t" << intercept << "\t" << lambda << "\t" << "\t" << sigma2 <<  "\t" << tgf << "\t" << alpha << "\t" <<ygf <<  endl;
                collect[0][j] = lambda;
                collect[1][j] = tgf;
                collect[2][j] = ygf;
                j++;
            }
        }	
//        if (i==5)
//            exit(1);
    }
    for (int i  = 0 ; i  < nef; i++) { 
            yl[i] = collect[1][i] * g[0];
            yu[i] = collect[1][i] * g[1];
    }
    if (io::debug >=1)
        cout << "nef =  " << nef << endl;
    gsl_sort (yl, 1, nef);
    gsl_sort (yu, 1, nef);
    int i1 = 0 ;
    int i2=  0 ;
    double cdf = 0 ;
    double c = 0 ;
    double bp = 0;
    double prevcdf = 0 ;
    double prevc = 0 ;
    double prevbp = 0;
    double w = (g[1]-g[0])*nef;
    double ygfl = -1;
    double ygfu = -1;
    while (1) { 
        prevc = c;
        prevcdf = cdf;
        prevbp = bp;
        if (i2 == nef)
            break;
        if (i1 < nef && yl[i1] < yu[i2]){
            c += 1/(yl[i1]/g[0]);
            bp = yl[i1];
            i1++;
        } else {
            c -= 1/(yu[i2]/g[1]);
            bp = yu[i2];
            i2++;
        }
//        cout << prevc << "\t" << cdf << "\t" << bp << "\t" << prevbp << endl;
        cdf += (prevc/w)*(bp-prevbp);
        if (cdf > 0.5 *tail &&  ygfl < 0 ) {
            ygfl = (0.5*tail - prevcdf )/(prevc/w) + prevbp;
        }
        if (cdf > 1 - 0.5*tail && ygfu < 0 ) { 
            ygfu = (1 - 0.5*tail - prevcdf )/(prevc/w) + prevbp;
        }
    }
    if (io::debug >=1 ) {
        cout << "cdf =  " << cdf << endl;
        cout << ygfl << "\t" << ygfu <<endl;
    }
    double mean, sd, l, u;
	statsfs.precision(0);
	statsfs.setf(ios::fixed);
    for (int  i =  0 ; i < ncollect; i ++) { 
        stats (collect[i], j, mean, sd, l , u);
        if (i==2)
            statsfs << noshowpoint <<  mean << "\t" << sd << "\t" << ygfl <<"\t" << ygfu <<endl;
        else
            statsfs << noshowpoint <<  mean << "\t" << sd << "\t" << l <<"\t" <<u <<endl;
    }
    statsfs.close ();
}

/*
double getposterior::findmle() {

    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;

    gsl_function F;

    F.function = &fit;
    F.params = this;

    double m = 1000, m_expected = 1200;
    double a = 0.0, b = 10000;
    double vala = fit (a,this);
    double valb = fit (b,this);
    double valm = fit (m,this);
    cout << vala << "\t" << valb <<"\t"<<valm << endl;

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);
    cout << "Here"<<endl;
    gsl_min_fminimizer_set (s, &F, m, a, b);
    cout << "Here"<<endl;

    printf ("using %s method\n",
            gsl_min_fminimizer_name (s));

    printf ("%5s [%9s, %9s] %9s %10s %9s\n",
            "iter", "lower", "upper", "min",
            "err", "err(est)");

    printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
            iter, a, b,
            m, m - m_expected, b - a);

    do
    {
        iter++;
        cout << "iter = " <<iter<<endl;
        status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        status 
            = gsl_min_test_interval (a, b, 0.001, 0.0);

        if (status == GSL_SUCCESS)
            printf ("Converged:\n");

        printf ("%5d [%.7f, %.7f] "
                "%.7f %+.7f %.7f\n",
                iter, a, b,
                m, m - m_expected, b - a);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_min_fminimizer_free (s);

    return status;
}*/

void getposterior::stats (double *x, int n, double &mean, double &sd, double &l, double &u) {
    gsl_sort (x, 1, n);
    mean = gsl_stats_mean (x,1,n);
    sd = sqrt(gsl_stats_variance (x,1,n));
    l = gsl_stats_quantile_from_sorted_data (x,1,n,0.5*tail);
    u = gsl_stats_quantile_from_sorted_data (x,1,n,1-0.5*tail);
}

int main (int argc, char* argv[]) {
    getposterior gp = getposterior (argc, argv);
    gp.findmle ();
    gp.mcmc ();
}
