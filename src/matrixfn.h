#ifndef MATRIXFN_H
#define MATRIXFN_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;
using namespace std;

namespace matrixfn  {

/* Matrix inversion routine.
 *  Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
    template<class T>
double invert(const matrix<T>& input, matrix<T>& inverse)
{
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<T> A(input);
    //
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());
    //
    // perform LU-factorization
    int res = lu_factorize(A, pm);
    double logdet = 0 ;
    for (int i  = 0 ; i < A.size1() ; i++) 
        logdet += log(abs(A(i,i)));



    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<T> (A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return logdet;
}
}
#endif

