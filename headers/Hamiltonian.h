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

#ifndef hamiltonian_h
#define hamiltonian_h
#include <stdio.h>
#include <complex>
#include <cmath>
#include <vector>
#include "print_matrices.h"

#include <iostream>
#include "./headers/create_annihilate.h"

#define PI 3.14159265359
extern string decay_type;
extern double coulomb_scaling;
extern double exp_scaling;
extern bool cluster;
extern bool pb;
extern int Nq;
extern int N;
extern double U;
extern double J;
extern vector<vector<double>> r;
extern bool EHM;
using namespace std;



double distance(vector<double> &r1,vector<double> &r2){
    double s = 0;
    for(int i = 0; i < r1.size(); ++i){
        s+=pow((r1[i] - r2[i]),2);
    }
    return sqrt(s);
}

double wavelet(double t,double t0,double sigma,double E,double wp){

    return E*exp(-pow((t-t0)/sigma,2)/2)*cos(wp*(t-t0));
    
}


int site_density(vector<vector<int>> &psi, int site){

    int temp = 0;
    if(psi[0][site] != 0){
        temp+=1;
    }
    if(psi[1][site] != 0){
        temp+=1;
    }
    return temp;
}


void h_non_eq(vector<complex<double>> &h,  double E,double wp,vector<vector<double>> &r, double t, double t0, double sigma){
    for(int i = 0; i < N; ++i){
        h[i*(N+1)] = (r[i][0])*wavelet(t,t0,sigma,E,wp);
//        h[i*(N+1)] = cos(PI*r[i][0]/2.0)*wavelet(t,t0,sigma,E,wp)
        }
}



void h_hopping(vector<complex<double>> &h,vector<vector<double>> &r,double J,double alpha,int N){
    double d;
    double dmax=distance(r[0],r[r.size()-1]);

    if(pb==false){
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                d=distance(r[i],r[j]);
                
                if(cluster == true){
                    if(i!=j){
                        h[i*N + j] = -J*exp(-alpha*(d - 1.0));
                    }
                }
                else if(cluster == false){
                    if(d <= 1.0 +1e-10 and d >= 1.0 - 1e-10){
                        h[i*N + j] = -J;
                    }
                }
            }
        }
    }
    else if(pb == true){
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                d=distance(r[i],r[j]);

                if(j != i and cluster == true){
                    
                    d = min(d,dmax+1.0-d);
                    h[i*N+j] = -J*exp(-alpha*(d-1.0));
                }
                else if(cluster == false){
                    if((d<= 1.0 + 1e-10 and d >= 1.0 - 1e-10) or (d <= N-1 + 1e-10 and d >= N-1 - 1e-10)){
                        h[i*N + j] = -J;
                    }
                }
            }
        }
    }

}

void h_sing_part(vector<complex<double>> &h,vector<complex<double>> &h_hopping,vector<vector<double>> &r,double offset,double t,double t0,double sigma,double E,double wp,double qe){
    vector<complex<double>> h_t(N*N);
    h_non_eq(h_t, E, wp,r, t, t0, sigma);
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            h[i*N + j] = h_hopping[i*N + j] + h_t[i*N+j];
            
            if(i==j){
                h[i*(N + 1)] += pow(-1,i)*offset;
                
            }
            if(i == j and i < Nq){
                
                h[i*(N + 1)] += qe;
                
            }
        
        }
    }
}

double repel_onsite(vector<vector<int>> &state,double U){
    vector<int> indi;
    vector<int> indj;
    for(int i = 0; i < state[0].size(); ++i){
        if(state[0][i] != 0){
            indi.push_back(i);
        }
        if(state[1][i] != 0){
            indj.push_back(i);
        }
    }
    double temp = 0;
    for(int i = 0; i < indi.size(); ++i){
        for(int j = 0; j < indj.size(); ++j){
            if(indi[i] == indj[j]){
                temp+=1;
            }
        }
    }
    return temp*U;
}


double repel_lng_rng(vector<vector<int>> &state,double U, vector<vector<double>> &r){
    vector<int> indi;
    vector<int> indj;
    double dmax = distance(r[0],r[r.size()-1]);
    double d;
    for(int i = 0; i < state[0].size(); ++i){
        if(state[0][i] != 0){
            indi.push_back(i);
        }
        if(state[1][i] != 0){
            indj.push_back(i);
        }
    }

    vector<int> ind(indi.size() + indj.size());
    for(int i = 0; i < indi.size(); ++i){
        ind[i] = indi[i];
        ind[i+indi.size()] = indj[i];

    }
    for(int i = 0; i < indj.size(); ++i){
        ind[i+indi.size()] = indj[i];

    }
    double temp = 0;
    for(int i = 0; i < ind.size(); ++i){
        for(int j = i; j < ind.size(); ++j){
            
            if(pb == true){
                d = distance(r[ind[i]],r[ind[j]]);
                d = min(d,dmax+1.0-d);
            }
            else if(pb == false){
                d = distance(r[ind[i]],r[ind[j]]);
            }
            if(d>1e-10){
                
                if(decay_type == "coulomb"){
                    temp+=U*coulomb_scaling/d;
                }
                else if(decay_type == "exp"){
                    temp+=U*exp(-exp_scaling*d);

                }
            }
        }
    }
    return temp;
}
void Ham(vector<complex<double>> &H, vector<complex<double>> &h_sp,vector<vector<vector<int>>> &states){
    fill(H.begin(),H.end(),0);
    int M = states.size();
    for(int i = 0; i < M;++i){
        for(int j = i; j < M; ++j){
            if(j==i){
                for(int k = 0; k < N; ++k){
                    int dens = site_density(states[i],k);
                    H[i*(M + 1)] += double(dens)*h_sp[k*(N + 1)];
                }
                
                H[i*(M + 1)] += repel_onsite(states[i],U);
                if(EHM == true){
                    H[i*(M + 1)] += repel_lng_rng(states[i],U,r);
                }
                
                
            }
            else{
                for(int k = 0; k < N; ++k){
                    for(int l = 0; l < N; ++l){
                        if(k!=l){
                            if(states[j][1] == states[i][1]){
                                double sign_up = 1;
                                vector<vector<int>> temp_up;
                                temp_up = states[i];
                                
                                sign_up *= annihilate(temp_up,k,0);
                                sign_up *= create(temp_up,l,0);
                                if(temp_up == states[j]){
                                    H[i*M + j] += sign_up*h_sp[k*N + l];
                                }
                            }
                            if(states[j][0] == states[i][0]){
                                
                                double sign_down = 1;
                                vector<vector<int>> temp_down;
                                temp_down = states[i];
                                
                                sign_down *= annihilate(temp_down,k,1);
                                sign_down *= create(temp_down,l,1);
                                if(temp_down == states[j]){
                                    H[i*M + j] += sign_down*h_sp[k*N + l];
                                }
                            }
                        }
                    }
                }
            }
            H[j*M + i] = conj(H[i*M + j]);
        }
    }
}


void update_Ham(vector<complex<double>> &H, vector<vector<double>> &r,vector<vector<vector<int>>> &states, double E,double wp, double t, double sigma, double t0, double dt){
    
    vector<complex<double>> h_t1(N*N);
    vector<complex<double>> h_t2(N*N);

    h_non_eq(h_t1, E, wp,r, t-dt, t0, sigma);
    h_non_eq(h_t2, E, wp,r, t, t0, sigma);
    int M = states.size();
    for(int i = 0; i < M;++i){
        for(int k = 0; k < N; ++k){
            int dens = site_density(states[i],k);
            H[i*(M + 1)] += (h_t2[k*(N+1)] - h_t1[k*(N+1)])*double(dens);
        }
            
    }
    
    
}



#endif /* Hamiltonian_h */
