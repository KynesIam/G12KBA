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

#ifndef pi_h
#define pi_h
#include "linear_algebra.h"
#include "print_matrix.h"
#include <vector>
#include <thread>
#include <complex>
#include "kron_delta.h"
using namespace std;
extern int Ns;
extern int Nb;

void Pixyxy(vector<complex<double>> &pi, vector<complex<double>> &rho_les,vector<complex<double>> &rho_grt, vector<complex<double>> &G2xxxx,vector<complex<double>> &G2xyxy, vector<complex<double>> &W,vector<complex<double>> &W_nonloc)
{
    int N = G2xxxx.size();
    vector<complex<double>> rho_d(N);
    vector<complex<double>> G2xxxxW(N);
    vector<complex<double>> WG2xxxx(N);
    vector<complex<double>> G2xyxyW_nonloc(N);
    vector<complex<double>> W_nonlocG2xyxy(N);
    vector<complex<double>> pi1(N);
    vector<complex<double>> pi2(N);
    vector<complex<double>> pi3(N);
    vector<complex<double>> pi4(N);
    for(int i = 0; i < N; ++i){
        rho_d[i] = rho_grt[i] - rho_les[i];
    }
    
#if NOTHREADS
#pragma omp parallel
    {
        matrix_mult(G2xxxxW,G2xxxx,W);
        matrix_mult(WG2xxxx,W,G2xxxx);
        matrix_mult(G2xyxyW_nonloc,G2xyxy,W_nonloc);
        matrix_mult(W_nonlocG2xyxy,W_nonloc,G2xyxy);
#pragma omp barrier
        matrix_mult(pi1,G2xxxxW,rho_d);
        matrix_mult(pi2,rho_d,WG2xxxx);
        matrix_mult(pi3,G2xyxyW_nonloc,rho_d);
        matrix_mult(pi4,rho_d,W_nonlocG2xyxy);
    }
#else
    thread th_G2xxxxW(matrix_mult,ref(G2xxxxW),ref(G2xxxx),ref(W));
    thread th_WG2xxxx(matrix_mult,ref(WG2xxxx),ref(W),ref(G2xxxx));
    thread th_G2xyxyW_nonloc(matrix_mult,ref(G2xyxyW_nonloc),ref(G2xyxy),ref(W_nonloc));
    thread th_W_nonlocG2xyxy(matrix_mult,ref(W_nonlocG2xyxy),ref(W_nonloc),ref(G2xyxy));
    th_G2xxxxW.join();
    th_WG2xxxx.join();
    th_W_nonlocG2xyxy.join();
    th_G2xyxyW_nonloc.join();
    
    thread th_pi1(matrix_mult,ref(pi1),ref(G2xxxxW),ref(rho_d));
    thread th_pi2(matrix_mult,ref(pi2),ref(rho_d),ref(WG2xxxx));
    thread th_pi3(matrix_mult,ref(pi3),ref(G2xyxyW_nonloc),ref(rho_d));
    thread th_pi4(matrix_mult,ref(pi4),ref(rho_d),ref(W_nonlocG2xyxy));
    
    th_pi1.join();
    th_pi2.join();
    th_pi3.join();
    th_pi4.join();
#endif
    
    for(int i = 0; i < N; ++i){
        pi[i] = -(pi2[i] + pi4[i] - pi1[i] - pi3[i]);
    }
}



void Pixyxy_sparse(vector<complex<double>> &pi, vector<complex<double>> &rho, vector<complex<double>> &G2xxxx,vector<complex<double>> &G2xyxy, vector<complex<double>> &W,vector<complex<double>> &W_nonloc)
{
    int N = G2xxxx.size();
    vector<complex<double>> rhoT(rho.size());
    vector<complex<double>> G2xxxxW(N);
    vector<complex<double>> WG2xxxx(N);
    vector<complex<double>> A1(N);
    vector<complex<double>> A2(N);
    vector<complex<double>> A3(N);
    vector<complex<double>> A4(N);
    
    vector<complex<double>> G2xyxyW_nonloc(N);
    vector<complex<double>> W_nonlocG2xyxy(N);
    vector<complex<double>> B1(N);
    vector<complex<double>> B2(N);
    vector<complex<double>> B3(N);
    vector<complex<double>> B4(N);
    
    for(int i = 0; i < rho.size(); ++i){
        rhoT[i] = conj(rho[i]);
    }
    
    
#if NOTHREADS
    matrix_mult(G2xxxxW,G2xxxx,W);
    matrix_mult(WG2xxxx,W,G2xxxx);
    matrix_mult(G2xyxyW_nonloc,G2xyxy,W_nonloc);
    matrix_mult(W_nonlocG2xyxy,W_nonloc,G2xyxy);
    
    sparse_mult1(A1,G2xxxxW,rhoT,'R');
    sparse_mult2(A2,G2xxxxW,rho,'R');
    sparse_mult1(A3,WG2xxxx,rhoT,'L');
    sparse_mult2(A4,WG2xxxx,rho,'L');
    
    sparse_mult1(B1,G2xyxyW_nonloc,rhoT,'R');
    sparse_mult2(B2,G2xyxyW_nonloc,rho,'R');
    sparse_mult1(B3,W_nonlocG2xyxy,rhoT,'L');
    sparse_mult2(B4,W_nonlocG2xyxy,rho,'L');
#else
    thread th_G2W(matrix_mult,ref(G2xxxxW),ref(G2xxxx),ref(W));
    thread th_WG2(matrix_mult,ref(WG2xxxx),ref(W),ref(G2xxxx));
    thread th_G2W_nonloc(matrix_mult,ref(G2xyxyW_nonloc),ref(G2xyxy),ref(W_nonloc));
    thread th_W_nonlocG2(matrix_mult,ref(W_nonlocG2xyxy),ref(W_nonloc),ref(G2xyxy));
    th_G2W.join();
    th_WG2.join();
    th_G2W_nonloc.join();
    th_W_nonlocG2.join();
    
    thread th_A1(sparse_mult1,ref(A1),ref(G2xxxxW),ref(rhoT),'R');
    thread th_A2(sparse_mult2,ref(A2),ref(G2xxxxW),ref(rho),'R');
    thread th_A3(sparse_mult1,ref(A3),ref(WG2xxxx),ref(rhoT),'L');
    thread th_A4(sparse_mult2,ref(A4),ref(WG2xxxx),ref(rho),'L');
    th_A1.join();
    th_A2.join();
    th_A3.join();
    th_A4.join();
    
    thread th_B1(sparse_mult1,ref(B1),ref(G2xyxyW_nonloc),ref(rhoT),'R');
    thread th_B2(sparse_mult2,ref(B2),ref(G2xyxyW_nonloc),ref(rho),'R');
    thread th_B3(sparse_mult1,ref(B3),ref(W_nonlocG2xyxy),ref(rhoT),'L');
    thread th_B4(sparse_mult2,ref(B4),ref(W_nonlocG2xyxy),ref(rho),'L');
    th_B1.join();
    th_B2.join();
    th_B3.join();
    th_B4.join();
#endif
    
    for(int i = 0; i < N; ++i){
        pi[i] = (A2[i] + B2[i] - A1[i] - B1[i])  - (A4[i] + B4[i] - A3[i] - B3[i]);
    }
}



#endif
