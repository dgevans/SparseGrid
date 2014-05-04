//
//  Chebyshev.h
//  SparseGrid
//
//  Created by David Evans on 1/3/14.
//  Copyright (c) 2014 David Evans. All rights reserved.
//

#ifndef __SparseGrid__Chebyshev__
#define __SparseGrid__Chebyshev__

#include <iostream>
#include <eigen3/Eigen/Dense>

double ChebyshevPolynomial(long n, double x);

Eigen::VectorXd ChebyshevPolynomialVector(long n, double x);

Eigen::MatrixXd ChebyshevPolynomialMatrix(long n, const Eigen::MatrixXd &X);

double ChebyshevRoots(long n, long j);


#endif /* defined(__SparseGrid__Chebyshev__) */
