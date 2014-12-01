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

MatrixXd overlap(std::vector<Shell>& shells){

    //Number of Orbitals is length of shells
    //Set up overlap matrix
    MatrixXd SS = MatrixXd::Zero(shells.size(), shells.size());

    //Populate overlap matrix.
    for(int i = 0; i < shells.size(); i++){
        for(int j = 0; j < shells.size(); j++){
            
            //Need to sum over inner contractions
            double s_tmp = 0.0;

            for(int p = 0; p < shells[i].contr.size(); p++){
                for( int q = 0; q < shells[j].contr.size(); q++){
                    
                    //cout << p << "\t" << shells[i].contr[p] << "\t" << q << "\t" << shells[j].contr[q] << endl;
                    //Guassian 2-center overlap integral is
                    //[pi/(alpha +
                    //beta)]^{3/2}Exp[-(alpha*beta)/(alpha+beta)|R_{A} -
                    //R_{B}|^{2}]
                    
                    double t_alpha = shells[i].alpha[p];
                    double t_beta = shells[j].alpha[q];
                    
                    double tmp1 = pow(M_PI / ( t_alpha + t_beta ), 3.0/2.0) ;
                    

                    double tmp2 = exp( (-1.*t_alpha*t_beta/(t_alpha + t_beta)) * (shells[i].origin - shells[j].origin).squaredNorm() );
                    
                    double norm_const1 = pow( 2*t_alpha/M_PI, 3./4.) ;
                    double norm_const2 = pow( 2*t_beta/M_PI, 3./4.);
                   
                    s_tmp += norm_const1*norm_const2*shells[i].contr[p]*shells[j].contr[q]*tmp1*tmp2;

                }
            }

            SS(i,j) = s_tmp;

        }
    }


    return SS;

}


MatrixXd kineticEnergy(std::vector<Shell>& shells){

    //Number of Orbitals is length of shells
    //Set up overlap matrix
    MatrixXd TT = MatrixXd::Zero(shells.size(), shells.size());

    //Populate kinetic energy matrix.
    for(int i = 0; i < shells.size(); i++){
        for(int j = 0; j < shells.size(); j++){
            
            //Need to sum over inner contractions
            double t_tmp = 0.0;

            for(int p = 0; p < shells[i].contr.size(); p++){
                for( int q = 0; q < shells[j].contr.size(); q++){
                    
                    double t_alpha = shells[i].alpha[p];
                    double t_beta = shells[j].alpha[q];
                   
                    double tmp = t_alpha*t_beta/(t_alpha + t_beta);
                    double tmp1 = pow(M_PI / ( t_alpha + t_beta ), 3.0/2.0) ;
                    
                    double tmp2 = exp( -1*tmp  * (shells[i].origin - shells[j].origin).squaredNorm() );
                    double tmp3 = 3 - 2*tmp*(shells[i].origin - shells[j].origin).squaredNorm();
                    double n_const_1 = pow( 2*t_alpha/M_PI, 3./4.) ;
                    double n_const_2 = pow( 2*t_beta/M_PI, 3./4.);
                   
                    t_tmp += n_const_1*n_const_2*tmp*tmp3*tmp1*tmp2*shells[i].contr[p]*shells[j].contr[q];

                }
            }

            TT(i,j) = t_tmp;
        }
    }

    return TT;
}

void nucElec(std::vector<Shell>& shells, std::vector<Atom>& atoms){

    //Compute with Fourier transform mathod
   
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
            double r2 = xij*xij + yij*yij + zij*zij;
            double r = sqrt(r2);
            enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
        }
    }
    cout << "\tNuclear repulsion energy = " << enuc << endl;

    //GENERATE BASIS
    std::vector<Shell> shells = make_sto3g_basis(atoms); 

    //COMPUTE OVERLAPS
    MatrixXd SS = overlap(shells);

    //Compute Kinetic Energy
    MatrixXd TT =kineticEnergy(shells);

    //Compute Nuclear-electron attraction
    nucElec(shells, atoms);





}



