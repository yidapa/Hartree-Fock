#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//DEBUGGING PACKAGES
#include <typeinfo>

//Eigen matrix algebra library Include Guard
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

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
    for (int i = 0; i < natom; i++){
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

        const double angstrom_to_bohr = 1 / 0.52917721092;
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




std::vector<Shell> make_sto3g_basis(std::vector<Atom>& atoms){

    //USE 3 primite Gaussians for STO-3G Look in Szabo and Ostlund for
    //coefficients and their overlaps to test if working
    Eigen::VectorXd alpha(3);
    alpha << 0.168856, 0.623913, 3.42525;
    Eigen::VectorXd contrcoef(3);;
    contrcoef << 0.444635, 0.535328, 0.154329;
    Eigen::VectorXd A(3);
    A << 0.0, 0.0, 0.0; //position of first Hydrogen
    Eigen::VectorXd B(3);
    B << 0.0, 0.0, 1.4; //position of second Hydrogen along z-axis

    
    std::vector<Shell> shells;

    shells.push_back( Shell{alpha, contrcoef, A } );
    shells.push_back( Shell{alpha, contrcoef, B } );


    return shells;
}



int main(int argc, char* argv[]){

    //const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    char* filename = argv[1];

    std::vector<Atom> atoms = read_geometry(filename);

    //Get Number of Electrons
    int nelectron = 0;
    for (int i = 0; i < atoms.size(); ++i){
        nelectron += atoms[i].atomic_number;
    }

    const int ndocc = nelectron / 2;

    //Now We can calculate Nuclear Energy
    double enuc = 0.0;
    for ( int i = 0; i < atoms.size(); i++){
        for (int j = i + 1; j < atoms.size(); j++){
            double xij = atoms[i].x - atoms[j].x;
            double yij = atoms[i].y - atoms[j].y;
            double zij = atoms[i].z - atoms[j].z;
            cout << xij << "\t" << yij << "\t" << zij << endl;
            double r2 = xij*xij + yij*yij + zij*zij;
            double r = sqrt(r2);
            enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
        }
    }
    cout << "\tNuclear repulsion energy = " << enuc << endl;

    //GENERATE BASIS
    std::vector<Shell> shells = make_sto3g_basis(atoms); 

    //COMPUTE OVERLAPS
    MatrixXd D = MatrixXd::Zero(2,2);

    cout << D << endl;





}



