//
//  linear_algebra.h
//  
//
//  Created by Cian Reeves on 4/22/23.
//

#ifndef linear_algebra_h
#define linear_algebra_h

#include <cblas.h>
#include "clapack.h"
//#include <mkl.h>//Uncomment if using MKL for diaganolization/matrix multiplication
//#include "mkl_lapacke.h" //Uncomment if using MKL for diaganolization/matrix multiplication
#include <stdio.h>
#include <vector>
#include <complex>
using namespace std;

void diagonalize(vector<complex<double>> &h1, vector<complex<double>> &evecs,vector<double> &evals)
{
    __CLPK_integer n=int(sqrt(h1.size()));
    __CLPK_integer LDA=n;
    __CLPK_integer lwork=-1;
    __CLPK_integer lrwork=-1;
    __CLPK_integer liwork=-1;

    __CLPK_integer info;
    vector<__CLPK_doublecomplex> a(h1.size());
    __CLPK_doublecomplex lwkopt;
    __CLPK_doublereal lrwkopt;
    __CLPK_integer liwkopt;



    for (int j=0; j<LDA; ++j){
        for(int i = 0; i<LDA;++i){
            a[i*LDA + j] = {double(h1[i*LDA+j].real()),double(h1[i*LDA + j].imag())};
        }
    }
    char jobz='V';
    char uplo='L';
    zheevd_(&jobz, &uplo,  &n,& *a.begin(), &LDA, & *evals.begin(),&lwkopt,&lwork,&lrwkopt,&lrwork,&liwkopt,&liwork,&info);

    lwork=(int)lwkopt.r;
    lrwork=(int)lrwkopt;
    liwork=(int)liwkopt;


    vector<__CLPK_doublecomplex> work(lwork);
    vector<__CLPK_doublereal> Rwork(lrwork);
    vector< __CLPK_integer> Iwork(max(1,liwork));

    info = zheevd_(&jobz, &uplo,  &n,& *a.begin(), &LDA, & *evals.begin(),& *work.begin(),&lwork,& *Rwork.begin(),&lrwork,& *Iwork.begin(),&liwork,&info);
    for (int j=0; j < LDA; ++j){
        for(int i = 0; i < LDA; ++i){
            evecs[i*n + j] = {a[i*LDA + j].r,a[i*LDA + j].i};
        }
    }
    if(info!=0)
        cout<<"Error in diagonalization\n";
}

//void diagonalize(vector<complex<double>> &A, vector<complex<double>> &evecs,vector<double> &evals)
//{
///* Function to perform double precision diagonalization of a Hermitian matrix using MKL.  Uncomment if using MKL
//    Args:
//        A: Vector storing matrix to be diagonalized
//        evecs: vector to store eigenvectors of A
//        evals: vector to store eigenvalues of A
// */
//    MKL_INT n = Ns*Nb, info;
//    vector<MKL_Complex16> a(A.size());
//
//    for (int j=0; j<n; ++j){
//        for(int i = j; i<n;++i){
//            a[i*n + j] = {A[i*n+j].real(),A[i*n + j].imag()};
//        }
//    }
//    info = LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'L', n,& *a.begin(), n, & *evals.begin());
//    for (int j=0; j < n; ++j){
//        for(int i = 0; i < n; ++i){
//            evecs[i*n + j] = {a[i*n+j].real,a[i*n + j].imag};
//        }
//    }
//
//    if(info!=0)
//        cout<<"Error in diagonalization\n";
//}

void matrix_mult(vector<complex<double>> &AB,vector<complex<double>> &A, vector<complex<double>> &B)
{
    char trans = 'N';
    complex<double> alpha,beta;
    alpha = 1.0;
    beta = 0.0;
    int N = int(sqrt(A.size()));
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *A.begin(),N,& *B.begin(),N,&beta,& *AB.begin(),N);
}


void matrix_vec_mult(vector<complex<double>> &AB,vector<complex<double>> &A, vector<complex<double>> &B)
{
    char trans = 'N';
    complex<double> alpha,beta;
    alpha = 1.0;
    beta = 0.0;
    int N = int(sqrt(A.size()));
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, N, 1,N, &alpha,& *A.begin(),N,& *B.begin(),1,&beta,& *AB.begin(),1);
}

complex<double> inner_product(vector<complex<double>> &A,vector<complex<double>> &B){
    if(A.size() != B.size()){
        cout<<"WARNING: Trying to compute inner product of vectors with different dimension.\n";
    }
    complex<double> ip = 0;
    for(int i = 0; i < A.size(); ++i){
        ip+=A[i]*conj(B[i]);
    }
    return ip;
}

complex<double> trace(vector<complex<double>> &G){
    complex<double> trG;
    int M=int(sqrt(G.size()));
    for(int k = 0; k < M; ++k){
        trG +=G[k*(M+1)];
    }
    return trG;
}
#endif /* linear_algebra_h */
