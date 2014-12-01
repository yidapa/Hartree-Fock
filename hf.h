/*****************************************
 *  Header file for hf.cc
 *  Contains Shell and Atom structures
 *  *************************************/

#ifndef EIGEN_DENSE_
#define EIGEN_DENSE_

#include <Eigen/Dense>

#endif /* Eigen/Dense */

#ifndef EIGEN_EIGENVALUES_
#define EIGEN_EIGENVALUES_

#include <Eigen/Eigenvalues>

#endif /* Eigen/Eigenvalues */



//using namespace Eigen;

struct Atom {
    int atomic_number;
    double x,y,z;
};

struct Shell {
    // contracted Gaussina = angular momentum + sph/cart flag + contracoeffs

    //struct Contraction{
    //    int l;
    //    bool pure;
    //    std::vector<LIBINT2_REALTYPE> coeff;
    //};
    //
    //
    Eigen::VectorXd alpha; //exponents
    Eigen::VectorXd contr; //contractions
    Eigen::VectorXd origin; //ORIGIN
};

//typedef Eigen::Matrix<double, Dynamic, Dynamic> MatrixXd;
 
