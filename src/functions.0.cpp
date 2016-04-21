#include "std.h"
#include "vectorfn.h"
#include "functions.h"

bool functions::fileExists(const std::string& fileName)
{
		std::fstream fin;
		fin.open(fileName.c_str(),std::ios::in);
		if( fin.is_open() )
		{
				fin.close();
				return true;
		}
		fin.close();
		return false;
}

double functions::l2distance (int n, double *x, double *y) {
	double d = 0 ;
	for (int i = 0; i < n ; i++){
		d += (x[i]-y[i])*(x[i]-y[i]);
	}
	d =  pow (d,0.5);
	return d;
}


// t - jackknifed estimates, m - weights, theta - original estimate
// returns jackknife corrected estimates and jackknife se
pair<double,double> functions::weightedjack (vector<double> &t, vector<double> &m, double theta) { 
    if (t.size () != m.size()){
        cerr << "In functions::weightedjack. Mismatch in lengths of t and m"<<endl;
        exit (1);
    }

    int g = 0  ;
    double n = vectorfn::sum(m);

    vector<double> res (m.size(),0);
    double a = 0;
    for (int i =  0 ; i < m.size(); i++) {
        if (m[i]<=0)
            continue;
        g++;
        double h = n/m[i];
        a += (1-1/h)*t[i];
        double r = theta * h - (h - 1) * t[i];
        res[i] = r;
    }

    if (g==0) {
        cerr << "In functions::weightedjack. Number of blocks ==0"<<endl;
        exit (1);
    }

    double tj = theta * g;
    tj -= a;
    double sj = 0 ;
    for (int i = 0 ; i < m.size(); i++) { 
        if (m[i]<=0)
            continue;
        double h = n/m[i];
        sj += pow((res[i]-tj),2)/(h-1);
    }
    sj /= g;
    sj = pow(sj,0.5);

    return  pair<double,double> (tj,sj);
}

void functions::tokenize(const string& str,
				vector<string>& tokens,
				const string& delimiters)
{
		// Skip delimiters at beginning.
		string::size_type lastPos = str.find_first_not_of(delimiters, 0);
		//         // Find first "non-delimiter".
		string::size_type pos     = str.find_first_of(delimiters, lastPos);

		while (string::npos != pos || string::npos != lastPos)
		{
				// Found a token, add it to the vector.
				tokens.push_back(str.substr(lastPos, pos - lastPos));
				// Skip delimiters.  Note the "not_of"
				lastPos = str.find_first_not_of(delimiters, pos);
				// Find next "non-delimiter"
				pos = str.find_first_of(delimiters, lastPos);
		}
}




void functions::sort (double *a, int length, int *index){
	for (int i = 0; i < length; i ++)
		index[i] = i;
	for (int i = 0;	i < length; i++){
		for (int j = 0; j < length ; j++){
			if ( a[i]  < a [j] ) {
				double tmp0  = a[i];
				a[i] = a[j];
				a[j] = tmp0;

				int tmp1 = index[i];
				index[i] = index[j];
				index[j] = tmp1;
			}
		}
	}
}


double functions::getsum (double** x, int l1, int l2){
	double max = getmax (x,l1,l2);
	double sum = 0;
	for (int i = 0; i < l1; i++){
		for (int j = 0; j < l2; j++){
			double tmp = x[i][j];
			if ( !std::isinf(tmp) && !std::isnan(tmp)){
				sum += exp (tmp - max);
			}
		}
	}
	
	sum = log (sum) + max;
	return sum;
}

double functions::getsum (double* x, int length){
	double max = getmax (x,length);
	double sum = 0;
	for (int i = 0; i < length; i++){
		double tmp = x[i];
		if ( !isinf(tmp) && !isnan(tmp)){
			sum += exp (tmp - max);
		}
	}

	sum = log (sum) + max;
	return sum;
}

void functions::normalize (double** x, int l1, int l2){
	double max = getmax (x,l1,l2);
	double sum = 0;
	for (int i = 0; i < l1; i++){
		for (int j = 0; j < l2; j++){
			double tmp = x[i][j];
			if ( !isinf(tmp) && !isnan(tmp)){
				sum += exp (tmp - max);
			}
		}
	}
	
	sum = log (sum) + max;

	for (int i = 0; i < l1; i++){
		for (int j = 0; j < l2; j++){
			double tmp = x[i][j];
			if ( !isinf(tmp) && !isnan(tmp)){
				x[i][j] = ( tmp - sum );
			} else {
				x[i][j] = 0;
			}
		}
	}
}

void functions::normalize (double* x, int length){
	double max = getmax (x,length);
	double sum = 0;
	for (int i = 0; i < length; i++){
		double tmp = x[i];
		if ( !isinf(tmp) && !isnan(tmp)){
			sum += exp (tmp - max);
		}
	}

	sum = log (sum) + max;

	for (int i = 0; i < length; i++){
		double tmp = x[i];
		if ( !isinf(tmp) && !isnan(tmp)){
			x[i] = ( tmp - sum );
		} else {
			x[i] = 0;
		}
	}

}

double functions::getmax (double* x, int length){
	double max;
	int flag = 0;
	for (int i = 0 ; i < length; i++){
		double tmp = x[i];
		if (!isinf(tmp) && !isnan (tmp)) {
			if (flag==0){
				max = tmp;
				flag = 1;
			} else if (tmp > max){
				max = tmp;
			}
		}
	}
	return max;
}

double functions::getmax (double** x, int l1, int l2){
	double max;
	int flag = 0;
	for (int i = 0 ; i < l1; i++){
		for (int j = 0 ; j < l2; j++){
			double tmp = x[i][j];
			if (!isinf(tmp) && !isnan (tmp)) {
				if (flag==0){
					max = tmp;
					flag = 1;
				} else if (tmp > max){
					max = tmp;
				}
			}
		}
	}
	return max;
}

pair<double,double> functions::meansd (double *x, int l) {
    double m2 = 0;
    double m1 = 0 ;
    for (int i = 0  ; i < l ; i++) {
        m1 += x[i];
        m2 += x[i]*x[i];
    }
    if (l>0) {
        m1 /= l ;
    }
    m2 = m2  -l* m1*m1;
    if (l>1) { 
        m2 /= (l-1);
    }
    return pair<double,double>(m1,m2);
}


bool functions::isnumeric( const char* pszInput, int nNumberBase )
{
    istringstream iss( pszInput );

    if ( nNumberBase == 10 )
    {
        double dTestSink;
        iss >> dTestSink;
    }
    else if ( nNumberBase == 8 || nNumberBase == 16 )
    {
        int nTestSink;
        iss >> ( ( nNumberBase == 8 ) ? oct : hex ) >> nTestSink;
    }
    else
        return false;

    // was any input successfully consumed/converted?
    if ( ! iss )
        return false;

    // was all the input successfully consumed/converted?
    return ( iss.rdbuf()->in_avail() == 0 );
}


