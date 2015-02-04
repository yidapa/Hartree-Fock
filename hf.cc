#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>

//DEBUGGING PACKAGES
#include <typeinfo>

//Eigen matrix algebra library Include Guard
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR> //For orthogonalization of basis

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
        else if (element_label == "He")
            Z = 2;
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


    //Use switch statement to generate the contraction coefficients.  Integrals
    //currently only handle s-orbitals

    std::vector<Shell> shells;
    Eigen::VectorXd alpha(3);
    Eigen::VectorXd contrcoef(3);
    Eigen::VectorXd A(3);
    Shell tmp;

    for(int zz = 0; zz < atoms.size(); zz++){

        switch (atoms[zz].atomic_number) {

            case(1):
                alpha << 0.16885540, 0.62391373, 3.42525091;
                contrcoef << 0.44463454, 0.53532814, 0.15432897;
                A << atoms[zz].x, atoms[zz].y, atoms[zz].z; //position of first Hydrogen
                tmp.alpha = alpha;
                tmp.contr = contrcoef;
                tmp.origin = A;
                shells.push_back( tmp );
                //cout << "\tLoaded H atom" << endl;
                break;
            case(2):
                //alpha << 0.31364979, 1.15892300, 6.36242139;
                alpha << 0.48084429026249986, 1.7766911481187495, 9.753934615874998;
                contrcoef << 0.44463454, 0.53532814, 0.15432897;
                A << atoms[zz].x, atoms[zz].y, atoms[zz].z; //position of first Hydrogen
                tmp.alpha = alpha;
                tmp.contr = contrcoef;
                tmp.origin = A;
                shells.push_back( tmp );
                //cout << "\tLoaded He atom" << endl;
                break;
            default:
                cout << "\tSorry! Can't handle atom with atomic number " << atoms[zz].atomic_number << endl;
 
        }

    }

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



MatrixXd nucElec(std::vector<Shell>& shells, std::vector<Atom>& atoms){


    MatrixXd Vnuc = MatrixXd::Zero(shells.size(), shells.size());

    //Compute with Fourier transform mathod
    for (int zz =0; zz < atoms.size(); zz++){ //loop over all atoms

        MatrixXd t_Vnuc = MatrixXd::Zero(shells.size(), shells.size());

        //compute repulsion matrix
        for(int i = 0; i < shells.size(); i++){
            for(int j = 0; j < shells.size(); j++){

                //Need to sum over inner contractions
                double t_tmp = 0.0;

                for(int p = 0; p < shells[i].contr.size(); p++){
                    for( int q = 0; q < shells[j].contr.size(); q++){

                        double t_alpha = shells[i].alpha[p];
                        double t_beta = shells[j].alpha[q];
 
                        double tmp1 = t_alpha*t_beta/(t_alpha + t_beta);
                        double gamma = t_alpha + t_beta;
                        double K = exp( -1*tmp1  * (shells[i].origin - shells[j].origin).squaredNorm() );

                        VectorXd Rp = (t_alpha*shells[i].origin + t_beta*shells[j].origin)/(t_alpha + t_beta);
                        
                        VectorXd Rc(3);
                        Rc << atoms[zz].x, atoms[zz].y , atoms[zz].z ; 

                        double n_const_1 = pow( 2*t_alpha/M_PI, 3./4.) ;
                        double n_const_2 = pow( 2*t_beta/M_PI, 3./4.);

                        double dist = (Rp - Rc).norm();
                        
                        if (dist == 0.0){
                            
                            t_tmp += -1*atoms[zz].atomic_number*n_const_1*n_const_2*shells[i].contr[p]*shells[j].contr[q]*2*K*M_PI/gamma;

                        } else{
                            t_tmp += -1*atoms[zz].atomic_number*n_const_1*n_const_2*shells[i].contr[p]*shells[j].contr[q]*K*pow(M_PI, 3./2)*erf(pow(gamma,1./2)*dist)/(pow(gamma,3./2)*dist);

                        }

                    }
                }

                t_Vnuc(i,j) = t_tmp;
            }
        }
        Vnuc += t_Vnuc;
    }

    return Vnuc;
}




MatrixXd TwoElecV(std::vector<Shell>& shells, std::vector<Atom>& atoms){

    //Set up two-electron repluslion matrix. goes as r^{4}
    MatrixXd VV = MatrixXd::Zero( pow(shells.size(),2), pow(shells.size(),2) );

    //Indexing is canoncial physics notation <kl|ij>
    
    for(int i = 0; i < shells.size(); i++){
        for(int j = 0; j < shells.size(); j ++){
            for(int k = 0; k < shells.size(); k++){
                for(int l = 0; l < shells.size(); l++){
                    //chem notation is [ki|lj]
                    //int,int rho{r1,p} rho{r2,q}/r12
                    //Vpq = int Vp(r2) rho{r2,Q} dr2
                    //express with Boys function.
                    
                    double t_VV = 0.0;

                    //Now for each we need a contraction over all gaussians
                    for(int p = 0; p < shells[i].alpha.size(); p++){
                        for(int q = 0; q < shells[j].alpha.size(); q++){
                            for(int r = 0; r < shells[k].alpha.size(); r++){
                                for( int s = 0; s < shells[l].alpha.size(); s++){
                                    
                                    double t_alpha = shells[i].alpha[p];
                                    double t_beta = shells[j].alpha[q];
                                    double t_gamma = shells[k].alpha[r];
                                    double t_delta = shells[l].alpha[s];

                                    double Kki = exp( -1.*t_alpha*t_gamma*(shells[k].origin - shells[i].origin).squaredNorm()/(t_alpha + t_gamma) );

                                    double Klj = exp( -1.*t_beta*t_delta*(shells[l].origin - shells[j].origin).squaredNorm()/(t_beta + t_delta) );

                                    VectorXd Rp = (t_alpha*shells[i].origin + t_gamma*shells[k].origin)/(t_alpha + t_gamma);
                                    VectorXd Rq = (t_beta*shells[j].origin + t_delta*shells[l].origin)/(t_beta + t_delta);
                                    
                                    double n__1 = pow( 2*t_alpha/M_PI, 3./4.) ;
                                    double n__2 = pow( 2*t_beta/M_PI, 3./4.);
                                    double n__3 = pow( 2*t_gamma/M_PI, 3./4.) ;
                                    double n__4 = pow( 2*t_delta/M_PI, 3./4.);
                                    
                                    double tmp = 0.0;
                                    //Evaluation of Boys function...
                                    if ( (Rp-Rq).norm() == 0){
                                        //F_{0} == 1
                                        tmp = n__1;
                                        tmp *= n__2*n__3*n__4;
                                        tmp *= shells[i].contr[p]*shells[j].contr[q]*shells[k].contr[r]*shells[l].contr[s];
                                        tmp *= 2*pow(M_PI, 5./2)*Kki*Klj/( (t_alpha+t_gamma)*(t_beta+t_delta)*pow(t_alpha+t_beta + t_gamma + t_delta,1./2) );
                                        
                                    } else {
                                        tmp = n__1;
                                        tmp *= n__2*n__3*n__4;
                                        tmp *= shells[i].contr[p]*shells[j].contr[q]*shells[k].contr[r]*shells[l].contr[s];

                                        tmp *= 2*pow(M_PI, 5./2)*Kki*Klj/( (t_alpha+t_gamma)*(t_beta+t_delta)*pow(t_alpha+t_beta + t_gamma + t_delta,1./2) );
                                        double T = (Rp - Rq).squaredNorm()*(t_alpha+t_gamma)*(t_beta+t_delta)/(t_alpha+t_beta+t_gamma + t_delta);
                                        tmp *= pow(M_PI/T,1./2)*erf(sqrt(T))/2.;
                                        


                                    }
                                    
                                    t_VV += tmp;
             
                                }
                            }
                        }
                    }
                
                    //cout << "<kl|ij>\t" << k+1 << l+1 << i+1 << j+1 << "\t" << t_VV << endl;
                    VV(i*shells.size()+ j, k*shells.size() + l) = t_VV;                
                    //if ( l == 1 && k == 1 && i == 1 && j == 1 ){
                    //    exit(0);
                    //}
                    //update matrix element

                }
            }
        }
    }

    return VV;
}


MatrixXd ConstructG(const Eigen::MatrixXd& D, const Eigen::MatrixXd& VVee){

    /***********************************
     *
     *          Construct G-matrix
     *
     *
     * ********************************/
   
    int L = D.rows();
    MatrixXd G = Eigen::MatrixXd::Zero(D.rows(), D.cols());
    for(int ii = 0; ii < D.rows(); ii++){ //row el1
        for(int jj = 0; jj < D.cols(); jj++){ //col el1 
            //printf("G_{%i}^{%i}\n",ii,jj);
            for(int rr = 0; rr < D.rows(); rr++){  //row el2
                for(int ss = 0; ss < D.cols(); ss++){ //col el2
                     
                    //G(ii,jj) += D(ss,rr)*(VVee(ii*L + rr, jj*L + ss) - 0.5*VVee(rr*L + ii, jj*L + ss) );
                    G(ii,jj) += D(rr,ss)*(VVee(ii*L + ss, jj*L + rr) - 0.5*VVee(ss*L + ii, jj*L + rr) );


                }
            }
        //std::cout << G << endl;
        }
    }

    return G;

}

int main(int argc, char* argv[]){

    //const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    char* filename = argv[1];

    std::vector<Atom> atoms = read_geometry(filename);

    const int bas_dim = atoms.size();

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
    //cout << "\tNuclear repulsion energy = " << enuc << endl;

    //GENERATE BASIS
    std::vector<Shell> shells = make_sto3g_basis(atoms); 

    //COMPUTE OVERLAPS
    MatrixXd SS = overlap(shells);
    
    cout << "\n\n\tOverlap Matrix" << endl;
     for(int i = 0; i < bas_dim; i++){
        for(int j = 0; j < bas_dim; j++){
            cout << "\t<" << j << "|" << i <<
                "> =\t" << SS(i,j) << endl;
        }
    }
    
    Eigen::IOFormat HeavyFmt(15,0,",",",\n", "[", "]", "[","]");
    cout << SS.format(HeavyFmt) << endl;


    //Compute Kinetic Energy
    MatrixXd TT =kineticEnergy(shells);

    
    cout << "\n\n\tKinetic Energy Matrix" << endl;
    for(int i = 0; i < bas_dim; i++){
        for(int j = i; j < bas_dim; j++){
            cout << "\t<" << j << "|" << i <<
                "> =\t" << TT(i,j) << endl;
        }
    }

    cout << TT.format(HeavyFmt) << endl;

    //Compute Nuclear-electron attraction
    MatrixXd VV = nucElec(shells, atoms);
    cout << "\n\n\tPotential Matrix" << endl;
    for(int i = 0; i < bas_dim; i++){
        for(int j = i; j < bas_dim; j++){
            cout << "\t<" << j << "|" << i <<
                "> =\t" << VV(i,j) << endl;
        }
    }
 
    cout << VV.format(HeavyFmt) << endl;



    //Core is sum of TT and VV
    MatrixXd Hcore = TT + VV;


    //COMPUTE 2-ELECTRON INTEGRALS
    MatrixXd VVee = TwoElecV(shells, atoms);
    cout << "\n" << VVee.format(HeavyFmt) << endl;

    

    int L = Hcore.rows();

    cout << "\n\n\tTwo Electron Integrals" << endl;
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            for(int k = 0; k < L; k++){
                for(int l = 0; l < L; l++){
         
            cout << "\t<" << k << "," << l << "|" << 
                i << "," << j << "> =\t" <<  
                VVee(i*L+j,k*L+l) << 
                "\t[" << k << "," <<  i << "|" << l << "," << 
                j << "] =\t" << VVee(i*L+j,k*L+l) << endl;  
                }
            }
        }
    }
        
    

    /**************************************************
     *              SCF Procedure
     *  1) generate orthogonalization matrix
     *  2) generate guess density.  Use Overlap matrix
     *  3) Build Fock Matrix
     *  4) transform Fock matrix
     *  5) diagonalize to get new orthonormal set.
     *  6) transform back to AO basis
     *  7) generate AO density
     *  8) check convergence go back to 3.
     **************************************************/

          
    ////Do symmetric orthogonalization by default 
    ////If eigenvalues below threshold Psi4 uses 1.0E-4 then use canonical
    ////orthogonalization (1.0E-7)
    ////
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigenss(SS);

    MatrixXd ssqrt = MatrixXd::Zero(eigenss.eigenvectors().rows(), eigenss.eigenvectors().cols() );

    for(int ii =0; ii < eigenss.eigenvectors().rows(); ii++)
        ssqrt(ii,ii) += pow(eigenss.eigenvalues()[ii], -0.5);

    MatrixXd X = eigenss.eigenvectors()*ssqrt;

    //
    ////Eigenvalues from solver sorted in increasing order Find all eigenvalues
    ////below 1.0E-7
   
    //int count = 0;
    //for(int ii =0; ii < eigenss.eigenvalues().size(); ii++)
    //    if (eigenss.eigenvalues()[ii] <= 1.0e-7) count += 1;
    //
    ////count is the number of eigenvalues less than cutoff. 
    //cout << count << endl;
    //
    //exit(0);


    //Generate Guess from core matrix

    MatrixXd F_prime = X.transpose()*Hcore*X;

    //Now solve for c' and epislon

    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(F_prime);

    MatrixXd C = X*eigensolver.eigenvectors();


    //Now form atomic orbital density matrix
    MatrixXd D = 2*C.leftCols(nelectron/2)*(C.leftCols(nelectron/2).transpose());

    cout << "\n\tInitial Density\n" << D << endl;


    MatrixXd G = ConstructG(D, VVee);


    MatrixXd Dold = D;


   
    int max_iter = 10000;
    int counter = 0;
    double energy = 0.0;
    while(counter < max_iter){

        Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(X.transpose()*(Hcore + G)*X);
        C = X*eigensolver.eigenvectors();
        D = 2*C.leftCols(nelectron/2)*C.leftCols(nelectron/2).transpose();
        
        //cout << "\t Density matrix norm " << (Dold - D).norm() << endl;

        //calculate energy
        for(int ii = 0; ii < Hcore.rows(); ii++){
            for(int jj = 0; jj < Hcore.cols(); jj++){
                energy += 0.5*D(ii,jj)*(2*Hcore(ii,jj) + G(ii,jj));
            }
        }

        //cout << "\tEnergy " << energy + enuc << endl;
        //cout << "\t P11 P12 P22\t" << D(0,0) << " " << D(0,1) << " " << D(1,1) << endl;

        if ( (Dold - D).norm() <= 1.e-8 ){
            cout << "\n\n\tSCF steps =\t" << counter << endl;
            cout << "\tNuclear Energy =\t" << enuc << endl;
            cout << "\tElectronic Energy =\t" << energy << endl;
            cout << "\tTotal Energy =\t" << energy + enuc << "\n" <<endl;
            break;
        }

        energy = 0.0;
        Dold = D;

        G = ConstructG(D, VVee);
        counter += 1;


    }

}
