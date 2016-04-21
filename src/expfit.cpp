#include "expfit.h"

expfit::expfit (int col , double l, double h, double scale )  { 
    this->col = col;
    this->l = l ;
    this->h= h;
    this->scale = scale;
    this->nmin = -1;
    this->nmax = -1;
}

void expfit::set_data ( vector<double> &tx, vector<double> &ty) { 

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

double neglogd (const gsl_vector *v, void *param) { 
    double val  = 0;
    int np = 0 ;
    expfit *gp = (expfit *)param;
    for( int i = gp->nmin ; i <= gp->nmax; i++){
        double tmp = gp->y[i]-gsl_vector_get(v,0)*exp(-gsl_vector_get(v,1)*gp->x[i]/gp->scale);
        val += tmp*tmp;
        np ++;
    }
    val = val *0.5/np;
    return (val);
}

double expfit::findmle(double &mleintercept, double &mlelambda, double &mlesigma2) {
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
}




void expfit::read_input (string filename) { 
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

/*
int main (int argc, char *argv[]) { 
}*/
