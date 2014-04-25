//
//  SparseGrid.cpp
//  SparseGrid
//
//  Created by David Evans on 1/4/14.
//  Copyright (c) 2014 David Evans. All rights reserved.
//

#include "SparseGrid.h"
#include "Chebyshev.h"
#include <iostream>

using namespace Eigen;


std::vector<VectorXl> interpolator::getSparseCombinations(long i_d)
{
    std::vector<VectorXl> ret;
    if (i_d == 0) {
        VectorXl temp = VectorXl::Zero(d);
        for (long i = 0; i <= mu; i++) {
            temp(0) = i;
            ret.push_back(temp);
        }
        return ret;
    }
    std::vector<VectorXl> ret_old = getSparseCombinations(i_d-1);
    for (long ib = 0; ib < ret_old.size(); ib++) {
        VectorXl temp = ret_old[ib];
        long temp_sum = temp.sum();
        for (long i = 0; i <= mu-temp_sum; i++) {
            temp(i_d) = i;
            ret.push_back(temp);
        }
    }
    return ret;
}

long m(long i)
{
    if(i < 0)
        return 0;
    if(i == 0)
        return 1;
    return pow(2,i)+1;
}

double evaluateBasisFunction(const VectorXl &bf, const MatrixXd &BM)
{
    double ret = 1;
    for (long i = 0; i < bf.rows(); i++) {
        ret *= BM(i,bf(i));
    }
    return ret;
}

void getBasisMatrix(const Eigen::VectorXd &x, Eigen::MatrixXd &BM)
{
    long n = BM.cols();
    for (long i = 0; i < BM.rows(); i++) {
        BM.row(i) = ChebyshevPolynomialVector(n-1, x(i)).transpose();
    }
}

VectorXd interpolator::SmolyakGrids(long i)
{
    long n = m(i)-m(i-1);
    VectorXd ret(n);
    if(i == 0)
        ret(0) = 0;
    else if(i == 1)
    {
        ret(0) = -1;
        ret(1) = 1;
    }else
    {
        for (long j = 0; j < n; j++) {
            ret(j) = ChebyshevRoots(m(i)-1, 2*j+1);
        }
    }
    return ret;
}

VectorXl interpolator::SmolyakFunctions(long i)
{
    VectorXl ret(m(i) - m(i-1));
    for (long j = m(i-1); j < m(i); j++) {
        ret(j-m(i-1)) = j;
    }
    return ret;
}

void interpolator::fillBasisFunctions(const std::vector<VectorXl> &SI)
{
    basis_functions.clear();
    std::vector<VectorXl> temp_bfs;
    for (long i =0; i < SI.size(); i++) {
        std::vector<VectorXl> x;
        for(long j =0; j < d; j++)
            x.push_back(SmolyakFunctions(SI[i](j)));
        temp_bfs = buildGrid(x);
        basis_functions.insert(basis_functions.end(), temp_bfs.begin(), temp_bfs.end());
    }
}

void interpolator::fillSmolyakGrid(const std::vector<Eigen::VectorXl> &SI)
{
    std::vector<VectorXd> Xvector;
    std::vector<VectorXd> temp_grid;
    for (long i =0; i < SI.size(); i++) {
        std::vector<VectorXd> x;
        for(long j =0; j < d; j++)
            x.push_back(SmolyakGrids(SI[i](j)));
        temp_grid = buildGrid(x);
        Xvector.insert(Xvector.end(), temp_grid.begin(), temp_grid.end());
    }
    
    Xbasis = MatrixXd(Xvector.size(),d);
    for (long ix = 0; ix < Xvector.size(); ix++) {
        Xbasis.row(ix) =Xvector[ix].transpose();
    }
}

interpolator::interpolator(long d_, long mu_):d(d_),mu(mu_)
{
    std::vector<VectorXl> SI = getSparseCombinations();
    fillBasisFunctions(SI);
    fillSmolyakGrid(SI);
    
    Mu = VectorXd::Zero(d);
    Sigma = MatrixXd::Identity(d, d);
}

interpolator::interpolator(long d_, long mu_, const VectorXd &Mu_, const MatrixXd &Sigma_):d(d_),mu(mu_),Mu(Mu_),Sigma(Sigma_)
{
    std::vector<VectorXl> SI = getSparseCombinations();
    fillBasisFunctions(SI);
    fillSmolyakGrid(SI);
    if ((Sigma.rows() != d) || (Sigma.cols() != d)) {
        throw "Sigma must be of dimension dxd";
    }
    if (Mu.rows() != d)
        throw "Mu must have d elemens";
}


IFunction interpolator::fit(const VectorXd &Y)
{
    long N = Xbasis.rows();
    MatrixXd Phi(N,N);
    MatrixXd BM(d,m(mu));
    for (long i = 0; i < N; i++) {
        getBasisMatrix(Xbasis.row(i), BM);
        for (long j = 0; j < N; j++) {
            Phi(i,j) = evaluateBasisFunction(basis_functions[j], BM);
        }
    }
    VectorXd c = Phi.fullPivLu().solve(Y);
    
    return IFunction(Mu, Sigma, c, basis_functions,mu);
}



IFunction::IFunction(const VectorXd &Mu_, const MatrixXd &Sigma_, const VectorXd &c_,const std::vector<VectorXl> &bfs_,long mu_): mu(mu_),Mu(Mu_),Sigma(Sigma_),basis_functions(bfs_),c(c_)
{
    d = Mu.rows();
    SigmaInv = Sigma.inverse();
}

VectorXd IFunction::eval(const rMatrixXd &X)
{
    MatrixXd TX;
    if (X.cols() == d)
    {
        TX= X.rowwise()-Mu.transpose();
    }else if (X.rows() == d)
    {
        TX = X.transpose().rowwise() - Mu.transpose();
    }else
    {
        throw "Cannot evaluate: X incorrect dimensions";
    }
    TX *= SigmaInv;
    
    MatrixXd Phi(TX.rows(),c.rows());
    
    MatrixXd BM(d,m(mu));
    for(long i = 0; i< Phi.rows(); i++)
    {
        getBasisMatrix(TX.row(i), BM);
        for (long j = 0; j < Phi.cols(); j++) {
            Phi(i,j) = evaluateBasisFunction(basis_functions[j],BM);
        }
    }
    return Phi*c;
}






