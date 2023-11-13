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

#ifndef ground_state_h
#define ground_state_h
#include <stdio.h>
#include <cmath>
#include "./linear_algebra.h"
#include <complex>
#include <vector>
using namespace std;
void prep_gs(vector<complex<double>> &psi, vector<complex<double>> &H){
    /*
        Finds ground state by diagonalizing the equilibrium Hamiltonian.
        Args:
            psi: vector to store the ground state wavefunction
            H: vector to store the equilibrium many-body Hamiltonian
    */
    cout<<"Preparing ground state\n";
    int M = int(sqrt(H.size()));
    vector<complex<double>> evecs(M*M);
    vector<double> evals(M);
    diagonalize(H, evecs,evals);
    if(evals[0] <= evals[1]+1e-10 and evals[0]>=evals[1]-1e-10){
        cout<<"WARNING: Degenerate ground state\n";
    }

    for(int i = 0; i < M; ++i){
        cout<<"Energy eigenvalue #"<<i+1<<":"<<evals[i]<<"\n";
    }
    for(int i = 0; i < int(sqrt(H.size())); ++i){
        
            psi[i] = evecs[i];
    }
}

#endif
