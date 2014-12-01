#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//DEBUGGING PACKAGES
#include <typeinfo>

//Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

//Libint Gaussian integrals library
#include <libint2.h>
#include <libint2/cxxapi.h>

//package includes
#include "hf.h"


using namespace Eigen;
using namespace std;


std::vector<Atom> read_geometry(char* filename) {
    
    std::cout << "Will read geometry from " << filename << std::endl;

    ifstream is(filename);
    assert(is.good());

    int natom;
    is >> natom;

    string str;
    getline(is, str);

    Atom atom;

    std::vector<Atom> atoms(natom);
    for (auto int i = 0; i < natom; i++){
        std::string element_label;
        double Z, x,y,z; //Element, xyz-cart
        is >> element_label >> x >> y >> z;
    
        if (element_label == "H")
            Z = 1;
        else if (element_label == "C")
          Z = 6;
        else if (element_label == "N")
          Z = 7;
        else if (element_label == "O")
          Z = 8;
        else if (element_label == "F")
          Z = 9;
        else if (element_label == "S")
          Z = 16;
        else if (element_label == "Cl")
          Z = 17;
        else {
          std::cerr << "read_dotxyz: element label \"" << element_label << "\" is not recognized" << std::endl;
          throw "Did not recognize element label in .xyz file";
        }
        
        atoms[i].atomic_number = Z;

        const auto double angstrom_to_bohr = 1 / 0.52917721092;
        atoms[i].x = x * angstrom_to_bohr;
        atoms[i].y = y * angstrom_to_bohr;
        atoms[i].z = z * angstrom_to_bohr;
    }

    //Check to see if we have read the file correctly
    //for(auto int i = 0; i < natom; i ++){
    //    cout << atoms[i].atomic_number << "\t" << "x- " << atoms[i].x << endl;
    //}

    return atoms;
}




void make_sto3g_basis(std::vector<Atom>& atoms){

    //no cxxapi.h as referenced by hartree-fock code.
    //must use Libint_t construct to keep information for one shell set of
    //integrals



    //Libint_t is the type of data structure used to pass basis function datat
    //to libint.  One Libint_t instance keeps information for one shell set of
    //integrals (or a vector of shell sets, if using vectoriztation)
    //TO evaluate integrals over contraction functions, allocate an array of
    //Libint_t objects, one for each combination of primitives
    //
    //

    //USE 3 primite Gaussians for STO-3G Look in Szabo and Ostlund for
    //coefficients and their overlaps to test if working
    double alpha[3] = {0.168856, 0.623913, 3.42525};
    double contrcoef[3] = {0.444635, 0.535328, 0.154329};
    double A[3] = {0.0, 0.0, 0.0}; //position of first Hydrogen
    double B[3] = {0.0, 0.0, 1.4}; //position of second Hydrogen along z-axis



}



int main(int argc, char* argv[]){

    //const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    char* filename = argv[1];

    std::vector<Atom> atoms = read_geometry(filename);

    //Get Number of Electrons
    auto int nelectron = 0;
    for (auto int i = 0; i < atoms.size(); ++i){
        nelectron += atoms[i].atomic_number;
    }
    const auto int ndocc = nelectron / 2;

    cout << "num electrons " << ndocc << endl;

    //Now We can calculate Nuclear Energy
    auto double enuc = 0.0;
    for ( auto int i = 0; i < atoms.size(); i++){
        for (auto int j = i + 1; j < atoms.size(); j++){
            auto double xij = atoms[i].x - atoms[j].x;
            auto double yij = atoms[i].y - atoms[j].y;
            auto double zij = atoms[i].z - atoms[j].z;
            cout << xij << "\t" << yij << "\t" << zij << endl;
            auto double r2 = xij*xij + yij*yij + zij*zij;
            auto double r = sqrt(r2);
            enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
        }
    }
    cout << "\tNuclear repulsion energy = " << enuc << endl;

    /***
     * Generate Basis
     ****/


    typedef unsigned int uint;
    // this initializes internal Libint data structures -- must happen once per
    // process
    LIBINT2_PREFIXED_NAME(libint2_static_init)();     
    //or use cxxapi 

    

    //maximum vector lenght of library.  set to 1 if no vectorization support
    //const unsigned int veclen = LIBINT2_MAX_VECLEN



    //libint library maximum angular momentum
    const unsigned int ammax = std::min(3, LIBINT2_MAX_AM_ERI);

    cout << ammax << endl;


    // initalize the libint integrals  library

    //COMPUTE OVERLAPS
    MatrixXd D = MatrixXd::Zero(2,2);

    cout << D << endl;





}



