//
//  Chebyshev.cpp
//  SparseGrid
//
//  Created by David Evans on 1/3/14.
//  Copyright (c) 2014 David Evans. All rights reserved.
//

#include "Chebyshev.h"
#include <math.h>

using namespace Eigen;

double ChebyshevPolynomial(long n, double x)
{
    //first check if n = 0 or 1
    if (n == 0)
        return 1.0;
    if (n == 1)
        return x;
    return 2*x*ChebyshevPolynomial(n-1, x) -
        ChebyshevPolynomial(n-2, x);
}
/*
 *Computes all then Chebyshev polynomials
 */
VectorXd ChebyshevPolynomialVector(long n, double x)
{
    //first check if n = 0 or 1
    VectorXd ret(n+1);
    ret(0) = 1;
    if (n == 0)
        return ret;
    ret(1) = x;
    if (n == 1)
        return ret;
    for (long i = 2; i <= n; i++) {
        ret(i) = 2*x*ret(i-1) - ret(i-2);
    }
    return ret;
}

MatrixXd ChebyshevPolynomialMatrix(long n, const MatrixXd &X)
{
    MatrixXd Y(X.rows(),X.cols());
    for (int i = 0; i< X.rows() ; i++) {
        for (int j = 0; j < X.cols(); j++)
            Y(i,j) = ChebyshevPolynomial(n, X(i,j));
    }
    return Y;
}

double ChebyshevRoots(long n, long j)
{
    if(n > 0)
        return -cos( M_PI * j/n);
    return 0;
}
