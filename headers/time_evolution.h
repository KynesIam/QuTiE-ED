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


#ifndef time_evolution_h
#define time_evolution_h
#include "./linear_algebra.h"
#include "./Hamiltonian.h"
#include "./observables.h"


void time_evolve_quench(vector<complex<double>> &hhop, vector<complex<double>> &psi, vector<vector<complex<double>>> &rhos,double q){
    /*
        Time evolves system in special case of a quench. Here the system is prepared in the ground state of on Hamiltonian and then
        time evolved with another one.  In this case the time evolution operator is only computed once
        Args:
            h_hopping: The kinetic energy portion of the Hamiltonian(single particle hopping terms)
            psi: Stores the wavefunction of the system
            rhos: stores the density matrix at each time step of the time evolution
    */
    vector<complex<double>> rho(N*N);
    vector<complex<double>> hsp(N*N);

    int M=states.size();
    vector<complex<double>> V(M*M);
    vector<complex<double>> H(M*M);

    h_sing_part(hsp,hhop,r,offset,0,t0,sigma,E,wp,q);
    Ham(H, hsp,states);
    evolution_operator(V,H,dt);

    for(int i = 0; i < n_steps; ++i){
        rho_t1t1(rho,psi,psi);
        rhos.push_back(rho);
        matrix_vec_mult(psi,V,psi);

        cout<<double(i)/(n_steps)*100<<"%\r";
        cout.flush();
    }
}

void time_evolve(vector<complex<double>> &hsp, vector<complex<double>> &psi, vector<vector<complex<double>>> &rhos,vector<vector<complex<double>>> &psis, vector<double> &Et){
    
    /*
        Time evolves system in general case of time dependent Hamiltonian. The system is prepared in the ground state of one Hamiltonian and then kicked by a
        time dependent field.  In this case the time evolution operator is only computed at every time step
        Args:
            
            psi: Stores the wavefunction of the system
            rhos: stores the density matrix at each time step of the time evolution
            rhos: stores the wavefunction at each time step of the time evolution
            energy: stores the energy at each time step of the time evolution

    */
    vector<complex<double>> rho(N*N);
    int M=states.size();
    vector<complex<double>> V(M*M);
    vector<complex<double>> H(M*M);
    
    Ham(H, hsp,states);
    for(int i = 0; i < n_steps; ++i){
        rho_t1t1(rho,psi,psi);
        rhos.push_back(rho);
        psis.push_back(psi);
        Et.push_back(energy(H,psi));

        
        update_Ham(H,r,states, E,wp, i*dt, sigma, t0,dt);
        evolution_operator(V,H,dt);
        matrix_vec_mult(psi,V,psi);
        cout<<double(i)/(n_steps)*100<<"%\r";
        cout.flush();
    }
}
#endif
