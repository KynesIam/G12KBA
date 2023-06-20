//
//
//   The routine(s) in this file are a part of the
//                     G12KBA
//   suite, developed 2022, and copyrighted
//   to the authors: Cian Reeves and Vojtech Vlcek
//   at the University of California, Santa Barbara
//   and Khaled Ibrahim
//   at Lawrence Berkeley National Lab, Berkeley.
//
//
//  If you use or modify any part of this routine
//  the header should be kept and unmodified.
//
//
//


#ifndef mat_mul_h
#define mat_mul_h
#include <vector>
#include <complex>
#include <iostream>
#include <thread>
#include <cblas.h>
#include "clapack.h"

extern int Ns;
extern int Nb;
extern int Ns2;
extern int Nb2;
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
void matrix_mult(vector<complex<double>> &AB,vector<complex<double>> &A, vector<complex<double>> &B)
{
    char trans = 'N';
    complex<double> alpha,beta;
    alpha = 1.0;
    beta = 0.0;
    int N = int(sqrt(A.size()));
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, N, N,N, &alpha,& *A.begin(),N,& *B.begin(),N,&beta,& *AB.begin(),N);
}


void get_subblock1(vector<complex<double>> &subblock, vector<complex<double>> &A,int Is, int Js,int Ib, int Jb)
{
    /*
     Function that finds the following NxN subblock of the matrix A.  S = A_{IJ}.....A_{IJ+N}
                                                                          .  .       .
                                                                          .     .    .
                                                                          .        . .
                                                                          A_{I+NJ}...A_{I+NJ+N}
     I,J represent Is,Ib,Js,Jb and N reperesnts Nb and Ns
     Args:
        subblock: Vector to store subblock described above
        A: Matrix from which subblock will be taken (M>Nb*Ns)
        Is,Js: Site index at which subblock starts
        Ib,Jb: Band index at which subblock starts
        Ns,Nb: Number of bands and sites in the model

     */

    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    int idx1 = IDX_4D(js,jb,is,ib);
                    int idx2 = IDX_8D(js,Js,jb,Jb,is,Is,ib,Ib);
                    subblock[idx1] = A[idx2];
                }
            }
        }
    }
}

void get_subblock2(vector<complex<double>> &subblock, vector<complex<double>> &A,int Is, int Js,int Ib, int Jb)
{
    /*
     Function that finds the following NxN subblock of the matrix A.  S = A_{IJ}A_{IJ+N}...A_{IJ+(N-1)N}
                                                                          .  .             .
                                                                          .        .       .
                                                                          .             .  .
                                                                          A_{I+N(N-1)J}....A_{I+N(N-1)J+(N-1)N}
     I,J represent Is,Ib,Js,Jb and N reperesnts Nb and Ns
     Args:
        subblock: Vector to store subblock described above
        A: Matrix from which subblock will be taken
        Is,Js: Site index at which subblock starts
        Ib,Jb: Band index at which subblock starts
        Ns,Nb: Number of bands and sites in the model

     */

    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    int idx1 = IDX_4D(js,jb,is,ib);
                    int idx2 = IDX_8D(Js,js,Jb,jb,Is,is,Ib,ib);
                    subblock[idx1] = A[idx2];
                }
            }
        }
    }
}

void sparse_mult1(vector<complex<double>> &result,vector<complex<double>> &A2, vector<complex<double>> &A1,char side)
{
    /*
     Function that calculates the matrix product (I_{NxN}kronA1_{NxN})(A2_{N^2xN^2}) or (A2_{N^2xN^2})(I_{NxN}kronA1_{NxN}) with O(N^5) operations. N ~ NsNb
     Args:
        result: Vector of length N^2xN^2 to hold matrix multiplication result
        A2: Vector of length N^2xN^2
        A1: Vector of length NxN
        side: char that determines whether to calculate (I_{NxN}kronA1_{NxN})(A2_{N^2xN^2}), side = 'L', or (A2_{N^2xN^2})(I_{NxN}kronA1_{NxN}), side = 'R'.
     */
    if(side == 'L'){
        
        vector<complex<double>> A1A2_sub(Ns*Nb*Ns*Nb);
        vector<complex<double>> subblock(Ns*Nb*Ns*Nb);

        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        
                        get_subblock1(subblock,A2,is,js,ib,jb);
                        
                        matrix_mult(A1A2_sub,A1,subblock);
                        
                        for(int ks = 0; ks < Ns; ++ks){
                            for(int kb = 0; kb < Nb; ++kb){
                                for(int ls = 0; ls < Ns; ++ls){
                                    for(int lb = 0; lb < Nb; ++lb){
                                        
                                        int idx1 = IDX_8D(ls,js,lb,jb,ks,is,kb,ib);
                                        int idx2 = IDX_4D(ls,lb,ks,kb);
                                        result[idx1] = A1A2_sub[idx2];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else if(side == 'R'){
        
        vector<complex<double>> A2_subA1(Ns*Nb*Ns*Nb);
        vector<complex<double>> subblock(Ns*Nb*Ns*Nb);

        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        
                        get_subblock1(subblock,A2,is,js,ib,jb);
                        matrix_mult(A2_subA1,subblock,A1);
                        
                        for(int ks = 0; ks < Ns; ++ks){
                            for(int kb = 0; kb < Nb; ++kb){
                                for(int ls = 0; ls < Ns; ++ls){
                                    for(int lb = 0; lb < Nb; ++lb){
                                        
                                        int idx1 = IDX_8D(ls,js,lb,jb,ks,is,kb,ib);
                                        int idx2 = IDX_4D(ls,lb,ks,kb);
                                        result[idx1] = A2_subA1[idx2];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else{cout<<"ERROR: Invalid choice in sparse matrix arg:side";}

}

void sparse_mult2(vector<complex<double>> &result,vector<complex<double>> &A2, vector<complex<double>> &A1,char side)
{
    
    /*
     Function that calculates the matrix product (A1_{NxN}kronI_{NxN})(A2_{N^2xN^2}) or (A2_{N^2xN^2})(A1_{NxN}kronI_{NxN}) with O(N^5) operations. N ~ NsNb
     Args:
        result: Vector of length N^2xN^2 to hold matrix multiplication result
        A2: Vector of length N^2xN^2
        A1: Vector of length NxN
        side: char that determines whether to calculate (A1_{NxN}kronI_{NxN})(A2_{N^2xN^2}), side = 'L', or (A2_{N^2xN^2})(A1_{NxN}kronI_{NxN}), side = 'R'.
     */
    if(side == 'L'){
        
        vector<complex<double>> A1A2_sub(Ns*Nb*Ns*Nb);
        vector<complex<double>> subblock(Ns*Nb*Ns*Nb);
        
        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        
                        get_subblock2(subblock,A2,is,js,ib,jb);
                        matrix_mult(A1A2_sub,A1,subblock);
                        
                        for(int ks = 0; ks < Ns; ++ks){
                            for(int kb = 0; kb < Nb; ++kb){
                                for(int ls = 0; ls < Ns; ++ls){
                                    for(int lb = 0; lb < Nb; ++lb){
                                        
                                        int idx1 = IDX_8D(js,ls,jb,lb,is,ks,ib,kb);
                                        int idx2 = IDX_4D(ls,lb,ks,kb);
                                        result[idx1] = A1A2_sub[idx2];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else if(side == 'R'){
        
        vector<complex<double>> A2_subA1(Ns*Nb*Ns*Nb);
        vector<complex<double>> subblock(Ns*Nb*Ns*Nb);
        
        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        
                        get_subblock2(subblock,A2,is,js,ib,jb);
                        matrix_mult(A2_subA1,subblock,A1);
                        
                        for(int ks = 0; ks < Ns; ++ks){
                            for(int kb = 0; kb < Nb; ++kb){
                                for(int ls = 0; ls < Ns; ++ls){
                                    for(int lb = 0; lb < Nb; ++lb){
                                        
                                        int idx1 = IDX_8D(js,ls,jb,lb,is,ks,ib,kb);
                                        int idx2 = IDX_4D(ls,lb,ks,kb);
                                        result[idx1] = A2_subA1[idx2];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    else{cout<<"ERROR: Invalid choice in sparse matrix arg:side";}

}

void matrix_subtract(vector<complex<double>> &result, vector<complex<double>> &A,vector<complex<double>> &B)
{
    if(A.size() != B.size()){cout<<"ERROR: Trying to Subtract Matrices of Incompatible Dimension";}
    for(int i = 0; i < A.size(); ++i){
        result[i] = A[i] - B[i];
    }
}

void commutator(vector<complex<double>> &comm,vector<complex<double>> &G,vector<complex<double>> &h)
{
        //Function to calculate the commutator between the HF hamiltonian and the 1 particle greens function.
    int N = G.size();
    vector<complex<double>> hG(N);
    vector<complex<double>> Gh(N);
#if    NOTHREADS
#pragma omp parallel
    matrix_mult(hG,h,G);
    matrix_mult(Gh,G,h);
#else
    thread th_hG(matrix_mult,ref(hG),ref(h),ref(G));
    thread th_Gh(matrix_mult,ref(Gh),ref(G),ref(h));
    th_hG.join();
    th_Gh.join();
#endif
    for(int i = 0; i < N; ++i){
        comm[i] = hG[i] - Gh[i];
    }
}

//void diagonalize(vector<complex<double>> &h1, vector<complex<double>> &evecs,vector<float> &evals)
//{
//    __CLPK_integer n=Ns*Nb;
//    __CLPK_integer LDA=Ns*Nb;
//    __CLPK_integer lwork=-1;
//    __CLPK_integer info;
//    vector<__CLPK_complex> a(h1.size());
//    __CLPK_complex wkopt;
//    vector<float> Rwork(max(1,3*n-2));
//
//    for (int j=0; j<LDA; ++j){
//        for(int i = 0; i<LDA;++i){
//            a[i*LDA + j] = {float(h1[i*LDA+j].real()),float(h1[i*LDA + j].imag())};
//        }
//    }
//    char jobz='V';
//    char uplo='L';
//    cheev_(&jobz, &uplo,  &n,& *a.begin(), &LDA, & *evals.begin(),&wkopt,&lwork,& *Rwork.begin(),&info);
//
//    lwork = (int)wkopt.r;
//    vector<__CLPK_complex> work(lwork);
//
//    info = cheev_(&jobz, &uplo,  &n,& *a.begin(), &LDA, & *evals.begin(),& *work.begin(),&lwork,& *Rwork.begin(),&info);
//    for (int j=0; j < LDA; ++j){
//        for(int i = 0; i < LDA; ++i){
//            evecs[i*n + j] = {a[i*LDA + j].r,a[i*LDA + j].i};
//        }
//    }
//    if(info!=0)
//        cout<<"Error in diagonalization\n";
//}

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



#endif
