#include "std.h"
#include "functions.h"
#include "io.h"
#include "matrixfn.h"

namespace ublas = boost::numeric::ublas;
        int debug = 0 ;
class mvn {
    public:
        mvn (string);
        double likelihood ( string file) ;
        int dim ;
        ublas::matrix<double> c;
        ublas::matrix<double> prec;
        ublas::vector<double> m;
        double lndet;
};

mvn::mvn (string filename) { 
	ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;
    dim = 0;

	while ( std::getline (inp, line)){
		linenum ++;
        std::vector<string> toks ;
        functions::tokenize(line, toks, " \t" );
        if (linenum==1)  {
            dim = toks.size ();
            if (debug >= 2) { 
                cout << "dim = " << dim <<endl;
                cout << line << endl;
            }
            m.resize(dim);
            c.resize(dim,dim);
        }

        std::vector<double> x (dim);
        for (int i = 0 ; i < dim; i++){ 
            x[i] = atof(toks[i].c_str());
        }

        for (int i = 0; i < dim; i++) {
            m(i) += x[i];
            for (int j = 0 ; j < dim; j++){ 
                c(i,j) += x[i]*x[j];
            }
        }

    }
    inp.close ();

    if (linenum > 0){ 
        for (int i = 0; i < dim; i++) 
            m(i) /= linenum;
        for (int i = 0; i < dim; i++) {
            for (int j = 0 ; j < dim; j++){
                c(i,j) /= linenum;
                c(i,j) -= m(i)*m(j);
            }
        }
    }
    if (debug >= 2)
        cout <<  " Covariance =  " << c << endl;

    prec.resize(dim,dim);
    lndet = matrixfn::invert (c, prec);

    if (debug >= 2) { 
        cout <<  " Precision  = " << prec << endl;
        cout << " lndet = " << lndet <<endl;
    }
}


double mvn::likelihood ( string filename) { 
    ifstream inp (filename.c_str());
	if (!inp.is_open()){
		cerr << "Error reading file "<< filename <<endl;
		exit(1);
	}
	string line;
	int linenum  = 0;

    ublas::vector<double> x(dim);
	while ( std::getline (inp, line)){
		linenum ++;
		io::println ("line = " + line, 2);
        std::vector<string> toks ;
        functions::tokenize(line, toks, " \t" );
        if (toks.size() != dim) { 
            cerr << "Badly formatted file " << filename << endl;
            cout << dim << "\t" << toks.size() << endl;
            exit(1);
        }
        for (int i = 0 ; i < dim; i++){ 
            x(i) = atof(toks[i].c_str());
        }
        double ll = -0.5 * inner_prod (x-m, prod(prec,x-m)) - lndet;
        cout << ll << endl;
    }
    inp.close ();
}

int main (int argc, char* argv[]) {
	if (argc > 1) { 
        mvn m = mvn (string(argv[1]));	
        if (argc > 2)  
            m.likelihood (string(argv[2]));
    } else {
        cerr << "Usage: " << argv[0] << "  file " <<endl;
    }
}
