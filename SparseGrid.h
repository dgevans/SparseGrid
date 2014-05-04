//
//  SparseGrid.h
//  SparseGrid
//
//  Created by David Evans on 1/4/14.
//  Copyright (c) 2014 David Evans. All rights reserved.
//

#ifndef __SparseGrid__SparseGrid__
#define __SparseGrid__SparseGrid__

#include <iostream>
#include <eigen3/Eigen/Core>
#include <vector>


namespace Eigen {
    typedef Matrix<long,Dynamic,1> VectorXl;
    typedef Matrix<long,Dynamic,Dynamic> MatrixXl;
    typedef Matrix<double,Dynamic,Dynamic,RowMajor> rMatrixXd;
}

class IFunction {
    //
    long mu;
    long d;
    
    //Hold the position of the interpolated grid
    Eigen::VectorXd Mu;
    Eigen::MatrixXd Sigma;
    Eigen::MatrixXd SigmaInv;
    
    //holds the list of combinations of basis functions
    std::vector<Eigen::VectorXl> basis_functions;
    
    //holds the vecor of coefficients of the basis function
    Eigen::VectorXd c;
    
    friend struct IFunction_pickle;
    
public:
    
    IFunction(){};
    
    IFunction(const Eigen::VectorXd &Mu, const Eigen::MatrixXd &Sigma, const Eigen::VectorXd &c,const std::vector<Eigen::VectorXl> &bfs, long mu);
    
    Eigen::VectorXd eval(const Eigen::rMatrixXd &X);
    
    Eigen::VectorXd evalMap( const Eigen::Map<Eigen::rMatrixXd> &X)
    {
        return eval(X);
    }
};

class interpolator
{
    long d; //dimension of the interpolator
    
    long mu; //degree of sparsity
    
    Eigen::VectorXd Mu;
    
    Eigen::MatrixXd Sigma;
    
    
    //holds the list of combinations of basis functions
    std::vector<Eigen::VectorXl> basis_functions;
    
    //holds the list of basis points
    Eigen::MatrixXd Xbasis;
    
    
    //get the combinations of i that form the sparse grid
    std::vector<Eigen::VectorXl> getSparseCombinations(){return getSparseCombinations(d-1);};
    std::vector<Eigen::VectorXl> getSparseCombinations(long i_d);
    
    //get the basis functions
    void fillBasisFunctions(const std::vector<Eigen::VectorXl> &SI);
    void fillSmolyakGrid(const std::vector<Eigen::VectorXl> &SI);
    
    //compute the smolyak gridpoint and the basis function for
    //a given integer i.
    Eigen::VectorXd SmolyakGrids(long i);
    Eigen::VectorXl SmolyakFunctions(long i);
    
    
public:
    
    interpolator(long d, long mu);
    
    interpolator(long d, long mu, const Eigen::VectorXd &Mu, const Eigen::MatrixXd &Sigma);
    
    Eigen::MatrixXd getXbasis(){return (Xbasis*Sigma).rowwise() + Mu.transpose();};
    
    IFunction fit(const Eigen::VectorXd &Y);
    

};



template <class vector_t>
std::vector<vector_t> buildGrid(const std::vector<vector_t> x) {
    long d = x.size();
    Eigen::VectorXl Ns(x.size());
    long N = 1;
    for (long i = 0; i < d; i++) {
        Ns(i) = N;
        N *= x[i].rows();
    }

    std::vector<vector_t> ret;
    for(long j = 0; j < N; j++)
    {
        vector_t temp(d);
        long jtemp = j;
        for(long i = d-1; i >=0; i--)
        {
            temp(i) = x[i](jtemp/Ns(i));
            jtemp = jtemp % Ns(i);
        }
        ret.push_back(temp);
    }
    return ret;
}
#endif /* defined(__SparseGrid__SparseGrid__) */
