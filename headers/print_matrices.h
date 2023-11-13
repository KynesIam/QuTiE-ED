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


#ifndef print_matrices_h
#define print_matrices_h
#include <stdio.h>
void print(vector<complex<double>> &H){
    for(int i = 0; i < int(sqrt(H.size())); ++i){
        for(int j = 0; j < int(sqrt(H.size())); ++j){
            cout<<H[i*int(sqrt(H.size())) + j]<<" ";
        }
        cout<<"\n";
    }
    cout<<"\n";

}
void print_basis(vector<vector<vector<int>>> &vec){
    for(int k = 0; k < vec.size(); ++k){
        cout<<"[[";
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < vec[0][0].size(); ++j){
                cout<<vec[k][i][j]<<",";
            }
            cout<<"] ";
        }
        cout<<"]";
        
        cout<<"\n";
    }

}

void print_vec(vector<complex<double>> &psi){
    for(int i = 0; i < psi.size(); ++i){
        cout<<psi[i]<<" ";
    }
    cout<<"\n";

}

void print_rintvec(vector<int> &psi){
    for(int i = 0; i < psi.size(); ++i){
        cout<<psi[i]<<" ";
    }
    cout<<"\n";

}
#endif /* print_matrices_h */
