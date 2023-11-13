//
//  write_to_file.h
//  
//
//  Created by Cian Reeves on 4/22/23.
//

#ifndef write_to_file_h
#define write_to_file_h

#include <fstream>
#include <vector>
#include <stdio.h>

using namespace std;
void write_to_rfile(string filename, vector<double> data)
{
    std::ofstream myfile;
    myfile.open(filename);
    for(int i = 0; i<(data.size()); i++){
            myfile << data[i] << '\n';
    }
    myfile.close();
}
void write_to_cfile(string filename, vector<complex<double>> data)
{
    std::ofstream myfile1;
    std::ofstream myfile2;

    myfile1.open("imag_" + filename);
    myfile2.open("real_" + filename);
    
    for(int i = 0; i<(data.size()); i++){
            myfile1 << data[i].imag() << '\n';
            myfile2 << data[i].real() << '\n';

    }
    myfile1.close();
    myfile2.close();

}

void write_to_rfile2D(string filename,vector<vector<double>> data)
{
    std::ofstream myfile1;
    myfile1.open(filename);

    for(int j = 0; j < data[0].size(); ++j){
        for(int i = 0; i<(data.size()); i++){
            myfile1 <<data[i][j]<<" ";
        }
        myfile1<<"\n";
    }
    myfile1.close();
}

void write_to_cfile2D(string filename,vector<vector<complex<double> > > data)
{
    std::ofstream myfile1;
    std::ofstream myfile2;
    myfile1.open("imag_" + filename);
    myfile2.open("real_" + filename);

    for(int j = 0; j < data[0].size(); ++j){
        for(int i = 0; i<(data.size()); i++){
            myfile1 <<data[i][j].imag()<<" ";
            myfile2 <<data[i][j].real()<<" ";
        }
        myfile1<<"\n";
        myfile2<<"\n";
    }
    myfile1.close();
    myfile2.close();
}

#endif /* write_to_file_h */
