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

#ifndef basis_states_h
#define basis_states_h
#include <cmath>
#include <vector>
using namespace std;

void swap(int& x, int& y)
{
    int temp = x;
    x = y;
    y = temp;
}


void permutations(vector<vector<int> >& res, vector<int> nums, int l, int h)
{
    if (l == h) {
        res.push_back(nums);
        return;
    }
 
    for (int i = l; i <= h; i++) {
        swap(nums[l], nums[i]);
        permutations(res, nums, l + 1, h);
        swap(nums[l], nums[i]);
    }
}

void unique_permutations(vector<vector<int>> &unique_perms, vector<int> &state){
    
    vector<int> cur_state;
    vector<vector<int>> perms;
    permutations(perms,state,0,state.size()-1);
    for(int i = 0; i < perms.size(); ++i){
        cur_state = perms[i];
        int temp=0;
        for(int j = 0; j < unique_perms.size(); ++j){
            if(cur_state == unique_perms[j]){
                temp+=1;
                break;
            }
        }
        if(temp == 0){
            unique_perms.push_back(cur_state);
        }
    }
}

void make_basis(vector<vector<vector<int>>> &states, vector<int> &state){
    /*
        
     */
    vector<vector<int>> unique_perms;
    unique_permutations(unique_perms,state);
    for(int i = 0; i < unique_perms.size(); ++i){
        for(int j = 0; j < unique_perms.size(); ++j){
            states.push_back({unique_perms[i],unique_perms[j]});
        }
    }
}

void make_asym_basis(vector<vector<vector<int>>> &states, vector<int> &state1,vector<int> &state2){
    vector<vector<int>> unique_perms1;
    vector<vector<int>> unique_perms2;

    unique_permutations(unique_perms1,state1);
    unique_permutations(unique_perms2,state2);

    for(int i = 0; i < unique_perms2.size(); ++i){
        for(int j = 0; j < unique_perms1.size(); ++j){
            states.push_back({unique_perms2[i],unique_perms1[j]});
        }
    }
}

#endif
