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


#ifndef greens_function_h
#define greens_function_h

#include <stdio.h>
#include "./create_annihilate.h"
#include "./linear_algebra.h"
#include "./basis_states.h"
#include "./Hamiltonian.h"
#include "./evolution_operator.h"
#include <vector>
#include <complex>
#include "./print_matrices.h"


extern vector<vector<vector<int>>> states;
extern vector<int> state;
extern vector<int> state_Np1;
extern vector<int> state_Nm1;

extern int N;
extern double offset;
extern double U;
extern double E;
extern double t;
extern double J;
extern double quench_strength;
extern double sigma;
extern double t0;
extern bool EHM;
extern double dt;
extern int dim;
extern vector<int> state;
extern int n_steps;
extern double alpha;
extern double wp;
extern vector<vector<double>> r;


void Glss_t1t2(vector<complex<double>> &G,vector<complex<double>> psi_t1,vector<complex<double>> psi_t2,vector<complex<double>> U_Nm1){
    /*
        Computes the lesser component of the Green's function at times t_1,t_2, ie G^<(t_1,t_2) = -i<psi(t_2)| c^\dagger U(t_2,t_1) c |psi(t_1)>
     
        Args:
            G: Vector to store the result of the Green's function
            psi_t1: The wavefunction at time t_1
            psi_t2: The wavefunction at time t_2
            U_Nm1: The portion of the time evolution operator acting on the N-1 particle subspace, where N is the total number of particles
     */
    vector<vector<vector<int>>> states_Nm1;
    complex<double> im = {0,1};
    make_asym_basis(states_Nm1,state,state_Nm1);//Makes a basis of states with N/2 particles in one spin subspace and and N/2-1 in the other(where N=total number)
    vector<complex<double>> psi_t1_temp(psi_t1.size());
    vector<complex<double>> psi_Nm1(states_Nm1.size());
    int M1 = states.size();
    int M2 = states_Nm1.size();
    int sign;
    
    vector<vector<int>> temp_N;
    vector<vector<int>> temp_Nm1;
    vector<vector<int>> zeros_N(states[0][0].size(),vector<int>(states[0][1].size()));
    vector<vector<int>> zeros_Nm1(states_Nm1[0][0].size(),vector<int>(states_Nm1[0][1].size()));

    for(int k = 0; k < N; ++k){
        for(int l = 0; l < N; ++l){
            psi_t1_temp = psi_t1;
            for(int i = 0; i < M1; ++i){//Loop creates the wavefunction after it is acted on by the annihilation operator
                temp_N = states[i];
                sign = annihilate(temp_N,k,0);
                for(int j = 0; j < M2; ++j){
                    if(temp_N == states_Nm1[j]){
                        psi_Nm1[j] = double(sign)*psi_t1_temp[i];
                        psi_t1_temp[i] = 0;
                        break;
                    }
                }
                if(psi_t1_temp[i] != 0.0){
                    psi_t1_temp[i] = 0.0;
                }
            }
            matrix_vec_mult(psi_Nm1,U_Nm1,psi_Nm1);
            for(int i = 0; i < M2; ++i){
                temp_Nm1 = states_Nm1[i];
                sign = create(temp_Nm1,l,0);
                for(int j = 0; j < M1; ++j){
                    if(temp_Nm1 == states[j]){
                        psi_t1_temp[j] = double(sign)*psi_Nm1[i];
                        psi_Nm1[i] = 0;
                        break;

                    }
                    
                }
                if(psi_Nm1[i] != 0.0){
                    psi_Nm1[i] = 0.0;
                }
            }
            G[l*N + k] = im*inner_product(psi_t2,psi_t1_temp);
        }
    }
}

void Ggrt_t1t2(vector<complex<double>> &G,vector<complex<double>> psi_t1,vector<complex<double>> psi_t2,vector<complex<double>> U_Np1){
    /*
        Computes the greater component of the Green's function at times t_1,t_2, ie G^>(t_1,t_2) = -i<psi(t_2)| c U(t_2,t_1) c^dagger |psi(t_1)>
     
        Args:
            G: Vector to store the result of the Green's function
            psi_t1: The wavefunction at time t_1
            psi_t2: The wavefunction at time t_2
            U_Np1: The portion of the time evolution operator acting on the N+1 particle subspace, where N is the total number of particles
     */
    
    vector<vector<vector<int>>> states_Np1;
    complex<double> im = {0,1};
    make_asym_basis(states_Np1,state,state_Np1);//Makes a basis of states with N/2 particles in one spin subspace and and N/2+1 in the other(where N=total number)
    vector<complex<double>> psi_t1_temp(psi_t1.size());
    vector<complex<double>> psi_Np1(states_Np1.size());
    int M1 = states.size();
    int M2 = states_Np1.size();
    int sign;
    
    vector<vector<int>> temp_N;
    vector<vector<int>> temp_Np1;
    vector<vector<int>> zeros_N(states[0][0].size(),vector<int>(states[0][1].size()));
    vector<vector<int>> zeros_Np1(states_Np1[0][0].size(),vector<int>(states_Np1[0][1].size()));

    for(int k = 0; k < N; ++k){
        for(int l = 0; l < N; ++l){
            psi_t1_temp = psi_t1;
            for(int i = 0; i < M1; ++i){
                temp_N = states[i];
                sign = create(temp_N,k,0);
                for(int j = 0; j < M2; ++j){
                    if(temp_N == states_Np1[j]){
                        psi_Np1[j] = double(sign)*psi_t1_temp[i];
                        psi_t1_temp[i] = 0;
                        break;

                    }
                    
                }
                if(psi_t1_temp[i] != 0.0){
                    psi_t1_temp[i] = 0.0;
                }
            }
            matrix_vec_mult(psi_Np1,U_Np1,psi_Np1);
            for(int i = 0; i < M2; ++i){
                temp_Np1 = states_Np1[i];
                sign = annihilate(temp_Np1,l,0);
                for(int j = 0; j < M1; ++j){
                    if(temp_Np1 == states[j]){
                        psi_t1_temp[j] = double(sign)*psi_Np1[i];
                        psi_Np1[i] = 0;
                        break;

                    }
                    
                }
                if(psi_Np1[i] != 0.0){
                    psi_Np1[i] = 0.0;
                }
            }
            G[l*N + k] = -im*inner_product(psi_t2,psi_t1_temp);
        }
    }
}

void rho_t1t1(vector<complex<double>> &G,vector<complex<double>> psi_t1,vector<complex<double>> psi_t2){
    /*
        Computes the density matrix at times t_1=t_2, ie rho^>(t_1=t_2) = <psi(t_2)| c U(t_2,t_1) c^dagger |psi(t_1)>
     
        Args:
            G: Vector to store the result of the Green's function
            psi_t1: The wavefunction at time t_1
            psi_t2: The wavefunction at time t_1
     */
    vector<vector<vector<int>>> states_Nm1;
    complex<double> im = {0,1};
    make_asym_basis(states_Nm1,state,state_Nm1);
    vector<complex<double>> psi_t1_temp(psi_t1.size());
    vector<complex<double>> psi_Nm1(states_Nm1.size());
    int M1 = states.size();
    int M2 = states_Nm1.size();
    int sign;
    vector<vector<int>> temp_N;
    vector<vector<int>> temp_Nm1;

    for(int k = 0; k < N; ++k){
        for(int l = 0; l < N; ++l){
            psi_t1_temp = psi_t1;
            for(int i = 0; i < M1; ++i){
                temp_N = states[i];
                sign = annihilate(temp_N,k,0);
                for(int j = 0; j < M2; ++j){

                    if(temp_N == states_Nm1[j]){
                        psi_Nm1[j] = double(sign)*psi_t1_temp[i];
                        psi_t1_temp[i] = 0;
                        break;
                    }
                }
                if(psi_t1_temp[i] != 0.0){
                    psi_t1_temp[i] = 0.0;
                }
            }
            
            for(int i = 0; i < M2; ++i){
                temp_Nm1 = states_Nm1[i];
                sign = create(temp_Nm1,l,0);
                for(int j = 0; j < M1; ++j){

                    if(temp_Nm1 == states[j]){
                        psi_t1_temp[j] = double(sign)*psi_Nm1[i];
                        psi_Nm1[i] = 0;
                        break;
                    }
                }
                if(psi_Nm1[i] != 0.0){
                    psi_Nm1[i] = 0.0;
                }

            }
            G[k*N + l] = inner_product(psi_t1_temp,psi_t2);
        }
    }
}


void Glss_full(vector<vector<complex<double>>> &G_full,vector<vector<complex<double>>> &psis){
    /*
        Scans through all t_1 and t_2 and computes G^< for all combinations
     
        Args:
            G_full: 2D vector to store G^< for all t_1 and t_2
            psis: Stores the wavefunction from each step of the time evolution
     */
    cout<<"Computing G^{less}(t_1,t_2) for all t_1<t_2\n";
    vector<complex<double>> G(N*N);
    vector<vector<vector<int>>> states_Nm1;
    make_asym_basis(states_Nm1,state,state_Nm1);

    int M = states_Nm1.size();
    vector<complex<double>> H_Nm1(M*M);
    vector<complex<double>> U_Nm1(M*M);
    vector<complex<double>> U_Nm1_temp(M*M);
    vector<complex<double>> eye(M*M);
    
    for(int i = 0; i < M; ++i){
        eye[i*(M+1)] = 1.0;
    }

    vector<complex<double>> hhop(N*N);
    vector<complex<double>> hsp(N*N);

    h_hopping(hhop,r,J,alpha,N);
    double count = 0;

    for(int i = 0; i < psis.size(); ++i){
        
        h_sing_part(hsp,hhop,r,offset, t,t0,sigma,E, wp,quench_strength);
        Ham(H_Nm1, hsp,states_Nm1);
        evolution_operator_dagger(U_Nm1,H_Nm1,dt);

        for(int j = i; j < psis.size(); ++j){
            if(i==j){
                
                Glss_t1t2(G,psis[i],psis[j],eye);
                G_full.push_back(G);
            }
            else{
                Glss_t1t2(G,psis[j],psis[i],U_Nm1);
                G_full.push_back(G);

                update_Ham(H_Nm1, r, states_Nm1,  E,wp, (j)*dt, sigma, t0, dt);
                evolution_operator_dagger(U_Nm1_temp,H_Nm1,dt);
                matrix_mult(U_Nm1,U_Nm1,U_Nm1_temp);
            }
            ++count;
            cout<<2.0*double(count)/(n_steps*(n_steps-1))*100<<"%\r";
            cout.flush();
        }
        
    }
}

void Ggrt_full(vector<vector<complex<double>>> &G_full,vector<vector<complex<double>>> &psis){
    /*
        Scans through all t_1 and t_2 and computes G^> for all combinations
     
        Args:
            G_full: 2D vector to store G^> for all t_1 and t_2
            psis: Stores the wavefunction from each step of the time evolution
    */
    cout<<"Computing G^{grt}(t_1,t_2) for all t_1<t_2\n";
    vector<complex<double>> G(N*N);
    vector<vector<vector<int>>> states_Np1;
    make_asym_basis(states_Np1,state,state_Np1);
    int M = states_Np1.size();
    vector<complex<double>> H_Np1(M*M);
    vector<complex<double>> U_Np1(M*M);
    vector<complex<double>> U_Np1_temp(M*M);
    vector<complex<double>> eye(M*M);
    
    for(int i = 0; i < M; ++i){
        eye[i*(M+1)] = 1.0;
    }

    vector<complex<double>> hhop(N*N);
    vector<complex<double>> hsp(N*N);
    h_hopping(hhop,r,J,alpha,N);
    double count = 0;

    for(int i = 0; i < psis.size(); ++i){
        h_sing_part(hsp,hhop,r,offset, t,t0,sigma,E, wp,quench_strength);
        Ham(H_Np1, hsp,states_Np1);


        for(int j = i; j < psis.size(); ++j){
            if(i==j){
                
                Ggrt_t1t2(G,psis[i],psis[j],eye);
                G_full.push_back(G);
            }
            else{
                
                Ggrt_t1t2(G,psis[i],psis[j],U_Np1);
                G_full.push_back(G);

                update_Ham(H_Np1,  r, states_Np1,  E,wp, (j)*dt, sigma, t0, dt);
                evolution_operator(U_Np1_temp,H_Np1,dt);
                matrix_mult(U_Np1,U_Np1,U_Np1_temp);
            }
            ++count;
            cout<<2.0*double(count)/(n_steps*(n_steps-1))*100<<"%\r";
            cout.flush();
        }
        
    }
}


void GR_full(vector<vector<complex<double>>> &GR,vector<vector<complex<double>>> &Glss_full,vector<vector<complex<double>>> &Ggrt_full){
    /*
     Computes retarded component of the Green's function from, G^R = G^> - G^<, for all t_1 and t_2
     Args:
        GR: 2D vector to store G^R for all t_1 and t_2

        Glss_full: 2D vector to store G^< for all t_1 and t_2
        Ggrt_full: 2D vector to store G^> for all t_1 and t_2

     */
    vector<complex<double>> GR_tmp(N*N);
    for(int i = 0; i < Glss_full.size(); ++i){
        fill(GR_tmp.begin(), GR_tmp.end(),0);
        for(int j = 0; j < Glss_full[0].size(); ++j){
            GR_tmp[j] = Ggrt_full[i][j] - Glss_full[i][j];
        }
        GR.push_back(GR_tmp);
    }
}

#endif
