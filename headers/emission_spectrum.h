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
#ifndef emission_spectrum_h
#define emission_spectrum_h
#include <vector>
#include <complex>
#include "./linear_algebra.h"
extern int N;
extern double dt;
extern int n_steps;
#define PI 3.14159265359
double gaussian(double t, double tp, double delta)
{
    return (1.0/(delta*sqrt(2.0*PI)))*exp(-pow((t-tp)/delta,2)/2);
}

void td_spectrum(vector<double> &I, vector<vector<complex<double>>> &G_t1t2, double delta, double tp, double w_max,double w_min)
{
    /*
        Computes the time dependent spectrum for a given time evolution described by G_t1t2(greater,lesser or retarded component).  Calculated as double integral over t1 and t2 of  S(t1-tp)S(t2-tp)e^{-iw(t1-t2)} Tr[G(t1,t2)]. where S is a gaussian
        Args:
            I: vector to store time-dependent spectrum
            G_t1t2: Complex vector to store Green's function at all points on t1, t2 grid
            delta: std-dev of Gaussian(S)
            tp: Represents time at which system is probed
            w_max: Maximum energy that spectrum is calculated up to
            w_min: Minimum energy that spectrum is calculated from
     */
    complex<double> im = {0,1};
    complex<double> trG=0;

    double w = w_min;
    double dw = .0125;

    while(w < w_max){
        
        cout<<"Computing emssion spectrum, percent complete:"<<double(w - w_min)/(w_max - w_min)*100<<"%\r";
        cout.flush();
        
        int idx=0;
        double t1 = 0;
        double t2 = 0;
        complex<double> I_temp=0;
        for(int i = 0; i < n_steps; ++i){
            t2=t1;
            for(int j = i; j < n_steps; ++j){
                trG=0;

                for(int k = 0; k < N; ++k){
                    trG += G_t1t2[idx][k*(N+1)];
                }
                if(i==j){
                    I_temp += dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(trG);
                }
                else{
                    I_temp += dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(trG*exp(im*(t1-t2)*w) - conj(trG)*exp(-im*(t1-t2)*w));
                }
                idx+=1;
                t2+=dt;
            }
            t1+=dt;
        }
        I.push_back(I_temp.imag());
        w += dw;
    }
}



void td_spectrum_scan(vector<vector<double>> &I_wtp, vector<vector<complex<double>>> &G_t1t2,  double delta, double w_max,double w_min)
{
    /*
        Scans through different probe times, tp, and creates a list containing date of energy vs probe time
     
        Args:
            I: vector to store time-dependent spectrum
            G_t1t2: Complex vector to store Green's function at all points on t1, t2 grid
            delta: std-dev of Gaussian(S)
            w_max: Maximum energy that spectrum is calculated up to
            w_min: Minimum energy that spectrum is calculated from
     */
    
    complex<double> im = {0,1};
    complex<double> trG=0;
    complex<double> trG2=0;
    double dw = .025;
    double tp = 0.0;
    double tp_max=n_steps*dt;
    double dtp = 5;
    while(tp<tp_max){
        vector<double> I;
        td_spectrum(I,G_t1t2, delta,tp,w_max,w_min);
        I_wtp.push_back(I);
        tp+=dtp;
        cout<<double(tp)/(tp_max)*100<<"%\r";
        cout.flush();
    }
}
#endif
