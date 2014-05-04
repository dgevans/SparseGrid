//
//  SparseGrid.cpp
//  SparseGrid
//
//  Created by David Evans on 1/3/14.
//  Copyright (c) 2014 David Evans. All rights reserved.
//
#include "NumpyConverter.hpp"
#include "Chebyshev.h"
#include "SparseGrid.h"
#include <boost/python.hpp>
#include <string>

namespace bp = boost::python;
using namespace Eigen;
using namespace std;

void translator(const char* x) {
    PyErr_SetString(PyExc_UserWarning, x);
}

struct IFunction_pickle : bp::pickle_suite
{
    
    static
    boost::python::tuple
    getinitargs(IFunction const& w)
    {
        return boost::python::make_tuple();
    }
    
    static
    boost::python::tuple
    getstate(IFunction const& F)
    {
        bp::list py_bf;
        for (long i =0; i < F.basis_functions.size(); i++)
        {
            py_bf.append(F.basis_functions[i]);
        }
        
        //finally return tuple of all thes objects
        return bp::make_tuple(F.Mu,F.Sigma,F.SigmaInv,py_bf,F.c,F.mu);
    }
    
    static
    void
    setstate(IFunction& F, boost::python::tuple state)
    {
        F.Mu = bp::extract<VectorXd>(state[0]);
        F.Sigma = bp::extract<MatrixXd>(state[1]);
        F.SigmaInv = bp::extract<MatrixXd>(state[2]);
        const bp::list &py_bf = bp::extract<bp::list>(state[3]);
        std::vector<VectorXl> &bfs = F.basis_functions;
        bfs.clear();
        for (long i = 0; i < bp::len(py_bf); i++)
            bfs.push_back(bp::extract<VectorXl>(py_bf[i]));
        F.c = bp::extract<VectorXd>(state[4]);
        F.mu = bp::extract<long>(state[5]);
        F.d = F.Mu.rows();
        /*interp.INFO = bp::extract<interpolator_INFO>(state[0]);
        interp.N = bp::extract<long>(state[1]);
        
        interp.bf.clear();
        bp::list bf_states = bp::extract<bp::list>(state[2]);
        for (int i =0; i < interp.N; i++) {
            bp::tuple bf_state = bp::extract<bp::tuple>(bf_states[i]);
            if (interp.INFO.types[i] == "spline")
                interp.bf.push_back(new basis_splines());
            else if(interp.INFO.types[i] == "hermite")
                interp.bf.push_back(new basis_hermite());
            else
                throw "in pickle set state: unkown type";
            interp.bf[i].load_state(bp::extract<VectorXd>(bf_state[0]), bp::extract<VectorXl>(bf_state[1]));
        }
        interp.c = bp::extract<VectorXd>(state[3]);*/
    }
};

BOOST_PYTHON_MODULE(SparseGrid)
{
    
    NumpyConverter::Register<MatrixXd>();
    
    NumpyConverter::Register<rMatrixXd>();
    
    NumpyConverter::Register<VectorXd>();
    
    NumpyConverter::Register<MatrixXl>();
    
    NumpyConverter::Register<VectorXl>();
    
    using namespace boost::python;
    
    register_exception_translator<
    char *>(&translator);
    
    def("T",ChebyshevPolynomial);
    def("T",ChebyshevPolynomialMatrix);
    def("Roots",ChebyshevRoots);
    
    class_<IFunction>("IFunction")
    .def_pickle(IFunction_pickle())
    .def("__call__",&IFunction::eval)
    .def("__call__",&IFunction::evalMap);
    
    class_<interpolator>("interpolator",init<long,long>())
    .def(init<long,long,const VectorXd&,const MatrixXd&>())
    .def("getX",&interpolator::getXbasis)
    .def("fit",&interpolator::fit);
}