//
//
//   The routine(s) in this file are a part of the
//                     QuTiE-ED
//   suite, developed 2023, and copyrighted
//   to the authors: Cian Reeves and Vojtech Vlcek
//   at the University of California, Santa Barbara
//
//
//  If you use or modify any part of this routine
//  the header should be kept and unmodified.
//
//
//

#ifndef dipole_h
#define dipole_h
#include <vector>
using namespace std;



void dipole(vector<vector<double>> &p,vector<vector<complex<double>>> &rho,vector<vector<double>> &r)
{
    /*
    Function to compute the dipole for a given set of coordinate vectors, for 3 dimensions the
    dipole is calculated along the z direction by default and for 1 dim along the one dimension available
    Args:
        p: vector to store the computed dipole, with length given by the total number of time steps
        rho: n_steps*(Ns*Ns*Nb*Nb) vector that holds the density matrix along the entire time evolution
        r: vector that stores coordinate vectors defining the lattice sites
    */
    int N = int(sqrt(rho[0].size()));
    for(int i = 0; i < rho.size(); ++i){
        for(int k = 0; k < N; ++k){
            for(int j = 0; j < 3; ++j)
            p[i][j] += r[k][j]*rho[i][k*(N+1)].real();

        }
    }
}


double energy(vector<complex<double>> &H,vector<complex<double>> &psi){
    vector<complex<double>> temp(psi.size());
    matrix_vec_mult(temp,H,psi);
    return inner_product(temp,psi).real();
}
#endif
