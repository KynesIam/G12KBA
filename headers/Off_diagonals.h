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

#ifndef Off_diagonals_h
#define Off_diagonals_h

#include "clapack.h"
//#include "mkl_lapacke.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>
#include "linear_algebra.h"
using namespace std;

extern int Ns;
extern int Nb;
extern int Ns2;
extern int Nb2;
extern double dt_fixed;
extern double t1,t2,t3;
extern int time_steps;
extern vector<vector<double>> r;
#define pi 3.14159265358979323846

void U_delta(vector<complex<double>> &result, vector<complex<double>> &h)
{
    complex<double> im={0,1};
    int N = int(sqrt(result.size()));
    double dt = dt_fixed;
    vector<complex<double>> evecs(N*N);
    vector<double> evals(N);

    diagonalize(h,evecs,evals);

    vector<complex<double>> temp(N*N);
    vector<complex<double>> diag_mat(N*N);
    vector<complex<double>> evecs_dag(N*N);

    for(int i = 0; i < N;++i){
        diag_mat[i*(N + 1)] = exp(-im*complex<double>(evals[i])*dt);
    }

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            evecs_dag[i*N + j] = conj(evecs[j*N + i]);
        }
    }

    matrix_mult(temp,evecs_dag,diag_mat);
    matrix_mult(result,temp,evecs);
}

void G_less_t1t2_full(vector<vector<complex<double>>> &G_lss_t1t2, vector<vector<complex<double>>> &G_grt_t1t2, vector<vector<complex<double>>> &rho_tot, vector<complex<double>> &W,vector<complex<double>> &hsp)
{
    double dt = dt_fixed;
    int N = rho_tot.size();
    int count = 0;
    complex<double> im={0,1};

    cout<<"Computing G lesser(t1,t2) at all points t1<t2\n";
    vector<complex<double>> h1(Ns*Ns*Nb*Nb);
    vector<complex<double>> U_op(Ns*Ns*Nb*Nb);
    vector<complex<double>> U_op_tmp(Ns*Ns*Nb*Nb);
    vector<complex<double>> temp_result(Ns*Ns*Nb*Nb);
    vector<complex<double>> temp_result2(Ns*Ns*Nb*Nb);

    vector<complex<double>> rho_tmp(Ns*Ns*Nb*Nb);
    vector<complex<double>> rho_grt_tmp(Ns*Ns*Nb*Nb);
    vector<complex<double>> im_eye(Ns*Ns*Nb*Nb);
    for(int i = 0; i < Ns*Nb;++i){
        im_eye[i*(Ns*Nb+1)] = im;
    }
    
    for(int i = 0; i < N; ++i){
        h_HF(h1,hsp,rho_tot[i],W,i*dt,Tp,wp,E,t0,r);
        for(int k = 0; k < rho_tot[i].size(); ++k){
            rho_tmp[k]=im*rho_tot[i][k];
            rho_grt_tmp[k]=im*rho_tot[i][k] - im_eye[k];
        }
        U_delta(U_op,h1);
        for(int j = (i); j < N; ++j){

            if(i == j){
                G_lss_t1t2.push_back(rho_tmp);
                G_grt_t1t2.push_back(rho_grt_tmp);
            }
            else{
                
                matrix_mult(temp_result,rho_tmp,U_op);
                matrix_mult(temp_result2,rho_grt_tmp,U_op);

                G_lss_t1t2.push_back(temp_result);
                G_grt_t1t2.push_back(temp_result2);

                h_HF(h1,hsp,rho_tot[j],W,j*dt,Tp,wp,E,t0,r);

                U_delta(U_op_tmp,h1);
                matrix_mult(U_op,U_op,U_op_tmp);
            }

            ++count;

            cout<<2.0*double(count)/(N*(N-1))*100<<"%\r";
            cout.flush();

        }
    }
}


//void integrate_h(vector<complex<double>> &result, vector<vector<complex<double>>> &rho_tot,vector<complex<double>> &W,int start_index, int end_index,double dt2)
//{
//    double dt=dt_fixed;
//    fill(result.begin(),result.end(),0);
//
//    if(end_index-start_index!=1){
//        for(int idx = start_index; idx < end_index; ++idx){
//
//            vector<complex<double>> h1(Ns*Ns*Nb*Nb);
//            h_HF_Hubbard(h1,rho_tot[idx],W,t1,t2,t3,0,epsilon,0,idx*dt,Tp,wp,0,t0);
//
//            if(idx==start_index){
//                for(int j = 0; j < result.size(); ++j){
//                    result[j]+=dt*h1[j]/2.0;
//                }
//            }
//            else if(idx == end_index-1){
//                for(int j = 0; j < result.size(); ++j){
//                    result[j]+=dt2*h1[j]/2.0;
//                }
//            }
//            else if(idx == end_index - 2){
//                for(int j = 0; j < result.size(); ++j){
//                    result[j]+=(dt2+dt)*h1[j]/2.0;
//                }
//            }
//            else{
//                for(int j = 0; j < result.size(); ++j){
//                    result[j]+=dt*h1[j];
//                }
//            }
//        }
//    }
//}
//
//void G_adv(vector<complex<double>> &result, vector<complex<double>> &h_int)
//{
//    complex<double> im={0,1};
//    int N = int(sqrt(result.size()));
//    double dt = dt_fixed;
//    vector<complex<double>> evecs(N*N);
//    vector<double> evals(N);
//
//    diagonalize(h_int,evecs,evals);
//
//    vector<complex<double>> temp(N*N);
//    vector<complex<double>> diag_mat(N*N);
//    vector<complex<double>> evecs_dag(N*N);
//
//    for(int i = 0; i < N;++i){
//        diag_mat[i*(N + 1)] = im*(exp(im*complex<double>(evals[i])));
//    }
//
//    for(int i = 0; i < N; ++i){
//        for(int j = 0; j < N; ++j){
//            evecs_dag[i*N + j] = conj(evecs[j*N + i]);
//        }
//    }
//
//    matrix_mult(temp,evecs_dag,diag_mat);
//    matrix_mult(result,temp,evecs);
//}
//

//void G_ret(vector<complex<double>> &result, vector<complex<double>> &h_int)
//{
//
//    complex<double> im={0,1};
//    int N = int(sqrt(result.size()));
//    vector<complex<double>> evecs(N*N);
//    vector<double> evals(N);
//    diagonalize(h_int,evecs,evals);
//
//    vector<complex<double>> temp(N*N);
//    vector<complex<double>> diag_mat(N*N);
//    vector<complex<double>> evecs_dag(N*N);
//
//    for(int i = 0; i < N;++i){
//        diag_mat[i*(N + 1)] = -im*exp(-im*complex<double>(evals[i]));
//    }
//
//    for(int i = 0; i < N; ++i){
//        for(int j = 0; j < N; ++j){
//            evecs_dag[i*N + j] = conj(evecs[j*N + i]);
//        }
//    }
//    matrix_mult(temp,evecs_dag,diag_mat);
//    matrix_mult(result,temp,evecs);
//}

//void G_t1t2(vector<vector<complex<double>>> &G_lss_t1t2,vector<vector<complex<double>>> &G_grt_t1t2,  vector<vector<complex<double>>> &rho_lss_tot, vector<complex<double>> &W,double dt2)
//{
//    double dt = dt_fixed;
//    int N = rho_lss_tot.size();
//
//    vector<complex<double>> temp_result1(Ns*Ns*Nb*Nb);
//    vector<complex<double>> temp_result2(Ns*Ns*Nb*Nb);
//
//    vector<complex<double>> rho_grt(Ns*Ns*Nb*Nb);
//
//    for(int time_index = 0; time_index < N; ++time_index){
//
//        vector<complex<double>> GA(Ns*Ns*Nb*Nb);
//        vector<complex<double>> h_int(Ns*Ns*Nb*Nb);
//
//
//        integrate_h(h_int,rho_lss_tot,W,time_index,N,dt);
//
//        G_adv(GA,h_int);
//
//
//        for(int j = 0; j < Ns; ++j){
//            for(int k = 0; k < Ns; ++k){
//                if(j == k){
//                    rho_grt[j*Ns + k] = rho_lss_tot[time_index][j*Ns + k] - 1.0;
//                }
//                else{
//                    rho_grt[j*Ns + k] = rho_lss_tot[time_index][j*Ns + k];
//                }
//            }
//        }
//
//        matrix_mult(temp_result1,rho_lss_tot[time_index],GA);
//        matrix_mult(temp_result2,rho_grt,GA);
//
//        G_lss_t1t2.push_back(temp_result1);
//        G_grt_t1t2.push_back(temp_result2);
//    }
//}


#endif
