#ifndef FUNCTIONS_H

#define FUNCTIONS_H

class functions{
	public:
	static bool fileExists (const std::string&);
 	static double l2distance (int n, double *x, double *y) ;

	static void tokenize(const string& str,	vector<string>& tokens,	const string& delimiters = " ");
	static void sort (double *a, int length, int *index);

	static void normalize (double* x, int length);
	static void normalize (double** x, int l1, int l2);
	static double getsum (double* x, int length);
	static double getsum (double** x, int l1, int l2);
	static double getmax (double** x, int l1, int l2);
	static double getmax (double* x, int length);
    static pair<double,double> meansd (double *x, int length);
    static void random_permutation (double *x, int n);
    static bool isnumeric( const char* pszInput, int nNumberBase = 10);
    static pair<double,double> weightedjack (vector<double> &x, vector<double> &w, double estimate);
};

#endif
