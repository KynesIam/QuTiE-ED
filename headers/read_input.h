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

#ifndef read_input_h
#define read_input_h
#include <stdio.h>
#include <vector>
using namespace std;
vector<vector<vector<int>>> states;

int N=4;
double offset=1.0;
double wp=0;
double U=1.0;
double E=2.5;
double t = 0;
double J = 1;
double sigma = .5;
double t0 = 6;
bool EHM=true;
bool cluster=false;
double dt = .02;
double exp_scaling=.7;
double coulomb_scaling=.5;
vector<int> state;
vector<int> state_Np1;
vector<int> state_Nm1;
bool pb=false;
int n_steps  = int(200.0/dt);
string decay_type="coulomb";
double alpha=.7;
vector<vector<double>> r;
int Nq=0;
double quench_strength=0;
string read_var_val(string filename,string var)
{
    /*
        Function to read in variable values from input file
        Args:
            filename: name of input file
            var: name of variable being read in
     */
    std::fstream myfile1;
    myfile1.open(filename);
    if(myfile1.fail()){
        cout<<"Missing input file 'input' \n";
    }
    int line_counter = 0;
    string line;
    while (getline(myfile1,line)){
        line_counter++;
        if (line.find(var) != string::npos){
            break;
        }
    }
    myfile1.close();
    
    std::fstream myfile2;
    myfile2.open(filename);
    string s;
    for(int i = 0; i < line_counter + 1; ++i){
        getline(myfile2,s);
    }
    if(s==""){
        cout<<"Missing value for" <<var<<"\n";
        return "0";
    }
    myfile2.close();
    return s;
}


void assign_vals()
{
    /*
         Function that initializes all variables that are read in from input file. Also
         initializes the vectors that determine the lattice geometry and the basis states
         that are used to compute the Hamiltonian
     */
    
    J=stod(read_var_val("input","hopping"));
    alpha=stod(read_var_val("input","hopping_decay"));

    U=stod(read_var_val("input","interaction_strength"));
    exp_scaling=stod(read_var_val("input","exp_scaling"));
    coulomb_scaling=stod(read_var_val("input","coulomb_scaling"));
    
    E=stod(read_var_val("input","pulse_strength"));
    t0=stod(read_var_val("input","pulse_midpoint"));
    wp=stod(read_var_val("input","pulse_frequency"));
    sigma=stod(read_var_val("input","Tp"));
    
    quench_strength=stod(read_var_val("input","quench_strength"));
    Nq=stoi(read_var_val("input","Nq"));

    N=stod(read_var_val("input","Ns"));
    
    
    offset=stod(read_var_val("input","offset"));

    n_steps=stoi(read_var_val("input","time_steps"));
    dt=stod(read_var_val("input","dt"));
    
    
    if(read_var_val("input","EHM") == "true" or read_var_val("input","EHM") == "True"){
        EHM=true;
    }
    else if(read_var_val("input","EHM") == "false" or read_var_val("input","EHM") == "False"){
        EHM=false;
    }
    else{
        cout<<"Invalid specification for EHM parameter. performing full calculation for onsite model.\n";
        EHM=false;
    }
    
    
    if(read_var_val("input","pb") == "true" or read_var_val("input","pb") == "True"){
        pb=true;
    }
    else if(read_var_val("input","pb") == "false" or read_var_val("input","pb") == "False"){
        pb=false;
    }
    else{
        cout<<"Invalid specification for pb parameter. performing calculation with open boundary conditions.\n";
        pb=false;
    }
    
    
    if(read_var_val("input","cluster") == "true" or read_var_val("input","cluster") == "True"){
        cluster=true;
    }
    else if(read_var_val("input","cluster") == "false" or read_var_val("input","cluster") == "False"){
        cluster=false;
    }
    else{
        cout<<"Invalid specification for cluster parameter. performing calculation for NN hopping.\n";
        cluster=false;
    }
    
    
        

//    if(read_var_val("input","q_type")=="full"){
//        q_type="full";
//    }
//    else if(read_var_val("input","q_type")=="pulse"){
//        q_type="pulse";
//    }
//    else if(read_var_val("input","q_type")=="none"){
//        q_type="none";
//        quench_strength=0;
//        Nq=0;
//    }
    
    
    if(read_var_val("input","decay_type")=="exp"){
        decay_type="exp";
    }
    else if(read_var_val("input","decay_type")=="coulomb"){
        decay_type="coulomb";
    }
    else{
        cout<<read_var_val("input","decay_type")<<"\n";
        cout<<"Invalid specification of decay type, defaulting to exponential decay\n";
        decay_type="exp";
    }

    for(int i = 0; i < int(N/2); ++i){
        state.push_back(1);
    }
    for(int i = 0; i < int(N/2); ++i){
        state.push_back(0);
    }
    for(int i = 0; i < int(N/2)+1; ++i){
        state_Np1.push_back(1);
    }
    for(int i = 0; i < int(N/2)-1; ++i){
        state_Np1.push_back(0);
    }
    for(int i = 0; i < int(N/2)-1; ++i){
        state_Nm1.push_back(1);
    }
    for(int i = 0; i < int(N/2)+1; ++i){
        state_Nm1.push_back(0);
    }
    for(int i = 0; i < N; ++i){
        r.push_back({(N/2) - i - .5,0,0});
    }
    
}


#endif /* read_input_h */
