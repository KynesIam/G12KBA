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

#ifndef psi_h
#define psi_h
#include "linear_algebra.h"
#include <vector>
#include <thread>
#include <complex>
using namespace std;

void Psi(vector<complex<double>> &psi, vector<complex<double>> &rho_les, vector<complex<double>> &rho_grt, vector<complex<double>> &W)
{
    int N = W.size();
    vector<complex<double>> Wrho_les(N);
    vector<complex<double>> Wrho_grt(N);
    vector<complex<double>> psi1(N);
    vector<complex<double>> psi2(N);
#if NOTHREADS
#pragma omp parallel
    {
        matrix_mult(Wrho_les,W,rho_les);
        matrix_mult(Wrho_grt,W,rho_grt);
#pragma omp barrier
        matrix_mult(psi1,rho_grt,Wrho_les);
        matrix_mult(psi2,rho_les,Wrho_grt);
    }
#else
    thread th_Wrho_les(matrix_mult,ref(Wrho_les),ref(W),ref(rho_les));
    thread th_Wrho_grt(matrix_mult,ref(Wrho_grt),ref(W),ref(rho_grt));
    th_Wrho_les.join();
    th_Wrho_grt.join();
    
    thread th_psi1(matrix_mult,ref(psi1),ref(rho_grt),ref(Wrho_les));
    thread th_psi2(matrix_mult,ref(psi2),ref(rho_les),ref(Wrho_grt));
    th_psi1.join();
    th_psi2.join();
#endif
    
    for(int i = 0; i < N; ++i){
        psi[i] = -(psi1[i] - psi2[i]);
    }
}

void Psi_X(vector<complex<double>> &psi, vector<complex<double>> &rho_les, vector<complex<double>> &rho_grt, vector<complex<double>> &W)
{
    int N = W.size();
    vector<complex<double>> Wrho_les(N);
    vector<complex<double>> WX(N);

    vector<complex<double>> Wrho_grt(N);
    vector<complex<double>> psi1(N);
    vector<complex<double>> psi2(N);
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    
    
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    for(int ks = 0; ks < Ns; ++ks){
                        for(int kb = 0; kb < Nb; ++kb){
                            for(int ls = 0; ls < Ns; ++ls){
                                for(int lb = 0; lb < Nb; ++lb){
                                    int idx1 = IDX_8D(js,is,jb,ib,ls,ks,lb,kb);
                                    int idx2 = IDX_8D(ls,is,lb,ib,js,ks,jb,kb);
                                    if(is==js && is==ks && is==ls && ib==jb && ib==kb && ib==lb){
                                        WX[idx1] = W[idx1];
                                        
                                    }
                                    else{
                                        WX[idx1] = (2.0*W[idx1] - W[idx2]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


#if NOTHREADS
#pragma omp parallel
    {
        matrix_mult(Wrho_les,WX,rho_les);
        matrix_mult(Wrho_grt,WX,rho_grt);
#pragma omp barrier
        matrix_mult(psi1,rho_grt,Wrho_les);
        matrix_mult(psi2,rho_les,Wrho_grt);
    }
#else
    thread th_Wrho_les(matrix_mult,ref(Wrho_les),ref(WX),ref(rho_les));
    thread th_Wrho_grt(matrix_mult,ref(Wrho_grt),ref(WX),ref(rho_grt));
    th_Wrho_les.join();
    th_Wrho_grt.join();
    
    thread th_psi1(matrix_mult,ref(psi1),ref(rho_grt),ref(Wrho_les));
    thread th_psi2(matrix_mult,ref(psi2),ref(rho_les),ref(Wrho_grt));
    th_psi1.join();
    th_psi2.join();
#endif
    
    for(int i = 0; i < N; ++i){
        if((psi1[i] - psi2[i]).real() > 0 and (psi1[i] - psi2[i]).imag() > 0){
            psi[i] = -(min((psi1[i] - psi2[i]).real(),.1),min((psi1[i] - psi2[i]).imag(),.1));
        }
        else if((psi1[i] - psi2[i]).real() > 0 and (psi1[i] - psi2[i]).imag() < 0){
            psi[i] = -(min((psi1[i] - psi2[i]).real(),.1),max((psi1[i] - psi2[i]).imag(),-.1));
        }
        if((psi1[i] - psi2[i]).real() < 0 and (psi1[i] - psi2[i]).imag() > 0){
            psi[i] = -(max((psi1[i] - psi2[i]).real(),-.1),min((psi1[i] - psi2[i]).imag(),.1));
        }
        if((psi1[i] - psi2[i]).real() < 0 and (psi1[i] - psi2[i]).imag() < 0){
            psi[i] = -(max((psi1[i] - psi2[i]).real(),-.1),max((psi1[i] - psi2[i]).imag(),-.1));
        }
        else{
            psi[i] = -(psi1[i] - psi2[i]);
        }
    }
}


void Psi_sparse(vector<complex<double>> &psi, vector<complex<double>> &rho, vector<complex<double>> &W)
{
    int N = W.size();
    vector<complex<double>> A1(N);
    vector<complex<double>> A2(N);
    vector<complex<double>> A3(N);
    vector<complex<double>> A4(N);
    vector<complex<double>> A5(N);
    vector<complex<double>> A6(N);
    vector<complex<double>> B1(N);
    vector<complex<double>> B2(N);
    vector<complex<double>> B3(N);
    vector<complex<double>> B4(N);
    vector<complex<double>> C1(N);
    vector<complex<double>> C2(N);
    vector<complex<double>> rhoT(rho.size());
                                    
    for(int i = 0; i < rho.size(); ++i){
        rhoT[i] = conj(rho[i]);
    }

#if NOTHREADS
#pragma omp parallel
    {
        sparse_mult1(A1,W,rhoT,'R');//WIrho_tmp
        sparse_mult1(A2,W,rhoT,'L');//Irho_tmpW
        sparse_mult2(A3,W,rho,'R');//Wrho_tmpI
        sparse_mult2(A4,W,rho,'L');//rho_tmpIW
#pragma omp barrier
        matrix_subtract(A5,A3,A1);//Wrho_tmpI - WIrho_tmp
        matrix_subtract(A6,A2,A4);//Irho_tmpW - rho_tmpIW
#pragma omp barrier
        sparse_mult1(B1,A6,rhoT,'R');
        sparse_mult1(B2,A3,rhoT,'L');
        sparse_mult1(B3,A5,rhoT,'L');
        sparse_mult2(B4,A1,rho,'L');
#pragma omp barrier
        sparse_mult2(C1,B1,rho,'R');
        sparse_mult2(C2,B3,rho,'L');
    }
#else
    
    thread th_A1(sparse_mult1,ref(A1),ref(W),ref(rhoT),'R');//WIrho_tmp
    thread th_A2(sparse_mult1,ref(A2),ref(W),ref(rhoT),'L');//Irho_tmpW
    thread th_A3(sparse_mult2,ref(A3),ref(W),ref(rho),'R');//Wrho_tmpI
    thread th_A4(sparse_mult2,ref(A4),ref(W),ref(rho),'L');//rho_tmpIW
    th_A1.join();
    th_A2.join();
    th_A3.join();
    th_A4.join();

    matrix_subtract(A5,A3,A1);//Wrho_tmpI - WIrho_tmp
    matrix_subtract(A6,A2,A4);//Irho_tmpW - rho_tmpIW

    thread th_B1(sparse_mult1,ref(B1),ref(A6),ref(rhoT),'R');
    thread th_B2(sparse_mult1,ref(B2),ref(A3),ref(rhoT),'L');
    thread th_B3(sparse_mult1,ref(B3),ref(A5),ref(rhoT),'L');
    thread th_B4(sparse_mult2,ref(B4),ref(A1),ref(rho),'L');
    th_B1.join();
    th_B2.join();
    th_B3.join();
    th_B4.join();

    thread th_C1(sparse_mult2,ref(C1),ref(B1),ref(rho),'R');
    thread th_C2(sparse_mult2,ref(C2),ref(B3),ref(rho),'L');
    th_C1.join();
    th_C2.join();
    
#endif
    
    for(int i = 0; i < N; ++i){
        psi[i] = (C1[i] + C2[i] + B4[i] - B2[i]);
    }
}


void Psi_X_sparse(vector<complex<double>> &psi, vector<complex<double>> &rho, vector<complex<double>> &W)
{
    int N = W.size();
    vector<complex<double>> WX(N);

    vector<complex<double>> A1(N);
    vector<complex<double>> A2(N);
    vector<complex<double>> A3(N);
    vector<complex<double>> A4(N);
    vector<complex<double>> A5(N);
    vector<complex<double>> A6(N);
    vector<complex<double>> B1(N);
    vector<complex<double>> B2(N);
    vector<complex<double>> B3(N);
    vector<complex<double>> B4(N);
    vector<complex<double>> C1(N);
    vector<complex<double>> C2(N);
    vector<complex<double>> rhoT(rho.size());
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    
    
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    for(int ks = 0; ks < Ns; ++ks){
                        for(int kb = 0; kb < Nb; ++kb){
                            for(int ls = 0; ls < Ns; ++ls){
                                for(int lb = 0; lb < Nb; ++lb){
                                    int idx1 = IDX_8D(js,is,jb,ib,ls,ks,lb,kb);
                                    int idx2 = IDX_8D(ls,is,lb,ib,js,ks,jb,kb);
                                    if(is==js && is==ks && is==ls && ib==jb && ib==kb && ib==lb){
                                        WX[idx1] = W[idx1];
                                        
                                    }
                                    else{
                                        WX[idx1] = (2.0*W[idx1] - W[idx2]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    for(int i = 0; i < rho.size(); ++i){
        rhoT[i] = conj(rho[i]);
    }

#if NOTHREADS
#pragma omp parallel
    {
        sparse_mult1(A1,WX,rhoT,'R');//WIrho_tmp
        sparse_mult1(A2,WX,rhoT,'L');//Irho_tmpW
        sparse_mult2(A3,WX,rho,'R');//Wrho_tmpI
        sparse_mult2(A4,WX,rho,'L');//rho_tmpIW
#pragma omp barrier
        matrix_subtract(A5,A3,A1);//Wrho_tmpI - WIrho_tmp
        matrix_subtract(A6,A2,A4);//Irho_tmpW - rho_tmpIW
#pragma omp barrier
        sparse_mult1(B1,A6,rhoT,'R');
        sparse_mult1(B2,A3,rhoT,'L');
        sparse_mult1(B3,A5,rhoT,'L');
        sparse_mult2(B4,A1,rho,'L');
#pragma omp barrier
        sparse_mult2(C1,B1,rho,'R');
        sparse_mult2(C2,B3,rho,'L');
    }
#else
    
    thread th_A1(sparse_mult1,ref(A1),ref(WX),ref(rhoT),'R');//WIrho_tmp
    thread th_A2(sparse_mult1,ref(A2),ref(WX),ref(rhoT),'L');//Irho_tmpW
    thread th_A3(sparse_mult2,ref(A3),ref(WX),ref(rho),'R');//Wrho_tmpI
    thread th_A4(sparse_mult2,ref(A4),ref(WX),ref(rho),'L');//rho_tmpIW
    th_A1.join();
    th_A2.join();
    th_A3.join();
    th_A4.join();

    matrix_subtract(A5,A3,A1);//Wrho_tmpI - WIrho_tmp
    matrix_subtract(A6,A2,A4);//Irho_tmpW - rho_tmpIW

    thread th_B1(sparse_mult1,ref(B1),ref(A6),ref(rhoT),'R');
    thread th_B2(sparse_mult1,ref(B2),ref(A3),ref(rhoT),'L');
    thread th_B3(sparse_mult1,ref(B3),ref(A5),ref(rhoT),'L');
    thread th_B4(sparse_mult2,ref(B4),ref(A1),ref(rho),'L');
    th_B1.join();
    th_B2.join();
    th_B3.join();
    th_B4.join();

    thread th_C1(sparse_mult2,ref(C1),ref(B1),ref(rho),'R');
    thread th_C2(sparse_mult2,ref(C2),ref(B3),ref(rho),'L');
    th_C1.join();
    th_C2.join();
    
#endif
    
    for(int i = 0; i < N; ++i){
        psi[i] = (C1[i] + C2[i] + B4[i] - B2[i]);
    }
    
}





#endif
