/*****************************************
 *  Header file for hf.cc
 *  Contains Shell and Atom structures
 *  *************************************/

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
    std::vector<double> alpha; //exponents
    //std::vector<Contraction> contr; //contractions
    //std::array<double, 3> O; //ORIGIN
};

//typedef Eigen::Matrix<double, Dynamic, Dynamic> MatrixXd;
 
