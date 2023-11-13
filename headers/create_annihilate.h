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


#ifndef create_annihilate_h
#define create_annihilate_h
#include <vector>
#include <cmath>
using namespace std;

double create(vector<vector<int>> &psi, int site, int spin){
    
    /*
        Updates psi by creating a particle at "site" with spin "spin" if one is not present and otherwise
        annihilates the state psi.  Updates state and returns the sign from fermionic anti-commutation relations
        Args:
            psi: state being acted on by creation operator
            site: site at which particle is being created
            spin: Determines which spin sector the creation operator acts on
     */
    int temp = 0;
    int N = psi[0].size();
    if(spin == 1){
        for(int k = 0; k < site; ++k){
            temp+=psi[spin][k];
        }
        for(int k = 0; k < N; ++k){
            temp+=psi[0][k];
        }
    }
    else if(spin == 0){
        for(int k = 0; k < site; ++k){
            temp+=psi[spin][k];
        }
    }
    
    if(psi[spin][site] == 0){
        psi[spin][site]=1;
    }
    else{
        for(int i = 0; i < psi.size(); ++i){
            for(int j = 0; j < psi[0].size(); ++j){
                psi[i][j] = 0;
            }
        }
    }

    return pow(-1,temp);

}

double annihilate(vector<vector<int>> &psi, int site, int spin){

    /*
        Updates psi by annihilating a particle at "site" with spin "spin" if one is present and otherwise
        annihilates the whole state psi.  Updates state and returns the sign from fermionic anti-commutation relations
        Args:
            psi: state being acted on by annihilation operator
            site: site at which particle is being annihilated
            spin: Determines which spin sector the annihilation operator acts on
     */
    
    
    int temp = 0;
    int N = psi[0].size();
    if(spin == 1){
        for(int k = 0; k < site; ++k){
            temp+=psi[spin][k];
        }
        for(int k = 0; k < N; ++k){
            temp+=psi[0][k];
        }
    }
    else if(spin == 0){
        for(int k = 0; k < site; ++k){
            temp+=psi[spin][k];
        }
    }
    
    if(psi[spin][site] == 1){
        psi[spin][site]=0;
    }
    else{
        for(int i = 0; i < psi.size(); ++i){
            for(int j = 0; j < psi[0].size(); ++j){
                psi[i][j] = 0;
            }
        }
    }

    return pow(-1,temp);
}
#endif
