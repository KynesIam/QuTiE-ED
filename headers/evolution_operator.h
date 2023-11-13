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


#ifndef evolution_operator_h
#define evolution_operator_h
#include <stdio.h>
#include <vector>
#include <complex>
#include <cmath>
#include "./linear_algebra.h"
using namespace std;

void evolution_operator(vector<complex<double>> &result, vector<complex<double>> &H,double dt)
{
    /*
        Computes the trotterized time evolution operator as the matrix exponential of the Hamiltonian by diagonalizing the Hamiltonian,
        exponentiating the eigenvalues and inverting the diagonalization
        Args:
            result: Vector to store the matrix exponential
            H: Vector to store the full Hamiltonian
            dt: The time step used in the time evolution
    */
    complex<double> im={0,1};
    int N = int(sqrt(result.size()));
    vector<complex<double>> evecs(N*N);
    vector<double> evals(N);
    
    diagonalize(H,evecs,evals);
    
    vector<complex<double>> temp(N*N);
    vector<complex<double>> diag_mat(N*N);
    vector<complex<double>> evecs_dag(N*N);
    
    for(int i = 0; i < N;++i){
        diag_mat[i*(N + 1)] = exp(-im*complex<double>(evals[i])*dt);
    }
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            evecs_dag[i*N + j] = conj(evecs[j*N + i]);
        }
    }

    matrix_mult(temp,evecs_dag,diag_mat);
    matrix_mult(result,temp,evecs);
}


void evolution_operator_dagger(vector<complex<double>> &result, vector<complex<double>> &H,double dt)
{
    /*
        Computes the Hermitian conjugate of the trotterized time evolution operator as the matrix exponential of the Hamiltonian by diagonalizing the Hamiltonian,
        exponentiating the eigenvalues and inverting the diagonalization
        Args:
            result: Vector to store the matrix exponential
            H: Vector to store the full Hamiltonian
            dt: The time step used in the time evolution
    */
    
    complex<double> im={0,1};
    int N = int(sqrt(result.size()));
    vector<complex<double>> evecs(N*N);
    vector<double> evals(N);
    
    diagonalize(H,evecs,evals);
    
    vector<complex<double>> temp(N*N);
    vector<complex<double>> diag_mat(N*N);
    vector<complex<double>> evecs_dag(N*N);
    
    for(int i = 0; i < N;++i){
        diag_mat[i*(N + 1)] = exp(im*complex<double>(evals[i])*dt);
    }
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            evecs_dag[i*N + j] = conj(evecs[j*N + i]);
        }
    }

    matrix_mult(temp,evecs_dag,diag_mat);
    matrix_mult(result,temp,evecs);
}
#endif
