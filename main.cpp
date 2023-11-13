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


#include <stdio.h>
#include <complex>
#include <cmath>
#include <vector>
#include <iostream>
#include "./headers/ground_state.h"
#include "./headers/greens_function.h"

#include "./headers/basis_states.h"
#include "./headers/Hamiltonian.h"
#include "./headers/write_to_file.h"
#include "./headers/evolution_operator.h"
#include "./headers/observables.h"
#include "./headers/read_input.h"
#include "./headers/print_matrices.h"
#include "./headers/emission_spectrum.h"
#include "./headers/time_evolution.h"
#include "./headers/basis_states.h"

#include <chrono>

using namespace std;


extern int N;
extern double offset;
extern double U;
extern double E;
extern double t;
extern double J;
extern double sigma;
extern double t0;
extern bool EHM;
extern double dt;
extern int dim;
extern vector<int> state;
extern int n_steps;
extern double alpha;
extern vector<int> state_Nm1;
extern vector<vector<double>> r;



int main(){
    
    auto start = chrono::steady_clock::now();
    assign_vals();
    make_basis(states,state);
    int M = states.size();
    
    vector<complex<double>> hhop(N*N);
    vector<complex<double>> hsp(N*N);
    vector<complex<double>> rho_gs(N*N);
    vector<complex<double>> psi(M);
    vector<complex<double>> H(M*M);
    
    vector<vector<complex<double>>> rhos;
    vector<vector<complex<double>>> Ggrt;
    vector<vector<complex<double>>> Glss;
    vector<vector<complex<double>>> GR;
    vector<vector<complex<double>>> psis;
    vector<double> I;
    vector<double> energy;
    vector<vector<double>> p(n_steps,vector<double>(3));

    
    h_hopping(hhop,r,J,alpha,N);
    h_sing_part(hsp,hhop,r,offset,0,t0,sigma,E,wp,0.0);
    
    cout<<"Single particle Hamiltonian at t = 0:\n";
    print(hsp);
    
    cout<<"Making Hamiltonian\n";
    Ham(H, hsp,states);
    
    prep_gs(psi,H);
    cout<<"\n";
    
    rho_t1t1(rho_gs,psi,psi);
    cout<<"Ground state density matrix:\n";
    print(rho_gs);
    
    
    cout<<"Beginning time stepping procedure\n";
    time_evolve(hsp,psi,rhos,psis,energy);
    cout<<"completed time stepping procedure\n";

    Glss_full(Glss,psis);
    Ggrt_full(Ggrt,psis);
    GR_full(GR,Glss,Ggrt);

    

    td_spectrum(I,GR, 15,n_steps*dt/2,6,-6.0);
    dipole(p,rhos,r);

    write_to_rfile2D("dipole.txt",p);
    write_to_cfile2D("density_matrix.txt",rhos);
    write_to_cfile2D("G_lss_t1t2.txt",Glss);
    write_to_cfile2D("G_grt_t1t2.txt",Ggrt);
    write_to_cfile2D("GR.txt",GR);
    write_to_rfile("emission.txt",I);
    write_to_rfile("energy.txt",energy);


    auto end_suc = chrono::steady_clock::now();
    cout<<"\n Execution Succesful, elapsed time:"<< chrono::duration<double>(end_suc - start).count()<< " sec \n";

    return 0;
}
