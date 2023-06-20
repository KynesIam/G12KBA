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

#ifndef rk_solver_h
#define rk_solver_h
#include "HF.h" //Holds various versions of  Hartree-Fock Hamiltonian
#include "psi.h"//Holds two particle occupation term
#include "pi.h" //Function for the polarization term in GW
#include "G1_G2_t_step.h"//Holds update steps for the one and two particle Green's function
#include "linear_algebra.h"
#include "collision_integral.h"//Collision integral, describes many body effects on one particle Green's function
#include "write_to_file.h"
#include "print_matrix.h"
#include "fermi_func.h"
#include "read_input.h"
#include "two_body.h"
#include "make_rho_lg.h"
#include "Self_Energy.h"
#include "Off_diagonals.h"
#include "twoPGF.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <thread>
#include <complex>
#include <chrono>
#include <string>

using namespace std;

extern complex<double> im;
extern bool HF;
extern bool EHM;
extern bool sparse;
extern double decay_rate;
extern double t1;
extern double t2;
extern double t3;
extern double P;
extern double U;
extern int time_steps;
extern double dt_fixed;
extern double quench_strength;
extern int Nq;
extern int Ns;
extern int Ns2;
extern int Nb;
extern int Nb2;
extern string q_type;
extern vector<double> epsilon;
extern double quench_rate;
extern double quench_on_midpoint;
extern double quench_off_midpoint;
extern double E;
extern double t0;
extern double Tp;
extern double wp;
extern string SE;
extern double G2_damping;




int RK4(vector<complex<double>> &rho,vector<complex<double>> &G2xxxx, vector<complex<double>> &G2xyxy,vector<vector<complex<double>>> &rho_diag,vector<complex<double>> &W,vector<complex<double>> &hsp,vector<vector<double>> &r)
{
    cout<<"Beginning Time Stepping Procedure\n";
    
    double tmax=double(time_steps)*dt_fixed;
    auto start = chrono::steady_clock::now();
    double RK_weight,dt,t,t_temp;
    double eps=1e-3;

    vector<complex<double>> I(Ns*Ns*Nb*Nb);

    vector<complex<double>> rho_les(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> rho_grt(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> h2(Ns2*Ns2*Nb2*Nb2);
    
    vector<complex<double>> W_non_loc(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> h2G2xxxx_comm(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> h2G2xyxy_comm(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> h1rho_comm(Ns*Ns*Nb*Nb);
 
    vector<complex<double>> psi_xxxx(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> psi_xyxy(Ns2*Ns2*Nb2*Nb2);
    
    vector<complex<double>> pi_xxxx(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> pi_xyxy(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> G2xxxx_tmp(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> G2xyxy_tmp(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> rho_tmp(Ns*Ns*Nb*Nb);
    vector<complex<double>> h(Ns*Ns*Nb*Nb);
    vector<double> evals_max;
    w_non_loc(W_non_loc,W);

    
    for(int t_index = 0; t_index < time_steps; ++t_index)
    {
        rho_diag.push_back(rho);

        vector<double> progress;
//        if(t_index%5000 == 0 and t_index != 0){
//            if(HF == true){
//                write_to_cfile2D("rho_tot_HF.txt",rho_diag);
//            }
//            if(HF == false){
//                write_to_cfile2D("rho_tot_GKBA.txt",rho_diag);
//            }
//        }
        vector<complex<double>> G2_full(G2xxxx.size());
        two_PGF(G2_full,G2xyxy,rho);
        vector<complex<double>> evecs(Ns2*Ns2*Nb2*Nb2);
        vector<double> evals(Ns2*Nb2);
        diagonalize(G2_full, evecs,evals);
        evals_max.push_back(evals[Ns2*Nb2-1]);
        for(int i = 0; i < Ns; ++i){
            double num_dens = 0;
            for(int a = 0; a < Nb; ++a){
                num_dens += rho[(i+a*Ns)*Nb*Ns + i + a*Ns].real();
            }
            
            if(num_dens> Nb + eps || isnan(num_dens) || num_dens < -eps){
                cout<<"ERROR:Unphysical Number Density on Site "<<(i+1)<<" Terminating Program. \nPercent Completed:"<<double(t)/(tmax)*100<<"% \n";
                if(HF == true){
                    write_to_cfile2D("rho_tot_HF.txt",rho_diag);
                }
                if(HF == false){
                    write_to_cfile2D("rho_tot_GKBA.txt",rho_diag);
                }
                write_to_rfile("G2_max_eval.txt",evals_max);

                return 1;
            }
        }
        
        vector<complex<double>> RK_result_rho(Ns*Ns*Nb*Nb);
        vector<complex<double>> RK_result_G2xxxx(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> RK_result_G2xyxy(Ns2*Ns2*Nb2*Nb2);
        
        vector<complex<double>> K1(Ns*Ns*Nb*Nb);
        vector<complex<double>> K2xxxx(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> K2xyxy(Ns2*Ns2*Nb2*Nb2);
        
        for(int k = 0; k<4; ++k){
            vector<complex<double>> h1(Ns*Ns*Nb*Nb);

            if(k == 1 or k == 2){
                dt = dt_fixed/2;
                RK_weight = 2;
            }
            else{
                dt = dt_fixed;
                RK_weight = 1;
            }
            
            if(k!=0){
                t_temp = t + dt;
            }
            else{
                t_temp = t;
            }

            for(int i = 0; i < Ns*Ns*Nb*Nb; ++i){
                rho_tmp[i] = rho[i] + K1[i]*dt;
            }
            
            for(int i = 0; i < Ns2*Ns2*Nb2*Nb2; ++i){
                G2xxxx_tmp[i] = G2xxxx[i] + K2xxxx[i]*dt;
                G2xyxy_tmp[i] = G2xyxy[i] + K2xyxy[i]*dt;
            }
        
            h_HF(h1,hsp,rho_tmp,W,t_temp,Tp,wp,E,t0,r);
            HF2(h2,h1);
            commutator(h1rho_comm,rho_tmp,h1);
            
#if NOTHREADS
            if(U!=0 && HF == false){

                commutator(h2G2xxxx_comm,G2xxxx_tmp,h2);
                commutator(h2G2xyxy_comm,G2xyxy_tmp,h2);
                collision_int(I,G2xxxx_tmp,G2xyxy_tmp,W);


                if(sparse==true){
                    if(SE=="scnd_brn"){
                        Psi_X_sparse(psi_xyxy,rho_tmp,W);
                        Psi_X_sparse(psi_xxxx,rho_tmp,W_non_loc);
                        
                    }
                    
                    if(SE=="GW"){
                        Psi_sparse(psi_xyxy,rho_tmp,W);
                        Psi_sparse(psi_xxxx,rho_tmp,W_non_loc);
                        Pixyxy_sparse(pi_xxxx,rho_tmp,G2xyxy_tmp,G2xxxx_tmp,W,W_non_loc);
                        Pixyxy_sparse(pi_xyxy,rho_tmp,G2xxxx_tmp,G2xyxy_tmp,W,W_non_loc);
                    }
                    
                }
                else if(sparse==false){
                    if(SE=="scnd_brn"){
                        make_rho_les(rho_les,rho_tmp,Ns,Nb);
                        make_rho_grt(rho_grt,rho_tmp,Ns,Nb);
                        
                        Psi_X(psi_xyxy,rho_les,rho_grt,W);
                        Psi_X(psi_xxxx,rho_les,rho_grt,W_non_loc);

                    }
                    else if(SE=="GW"){
                        make_rho_les(rho_les,rho_tmp,Ns,Nb);
                        make_rho_grt(rho_grt,rho_tmp,Ns,Nb);
                        
                        Psi(psi_xyxy,rho_les,rho_grt,W);
                        Psi(psi_xxxxrho_les,rho_grt,W_non_loc);
                        Pixyxy(pi_xxxx,rho_les,rho_grt,G2xyxy_tmp,G2xxxx_tmp,W,W_non_loc);
                        Pixyxy(pi_xyxy,rho_les,rho_grt,G2xxxx_tmp,G2xyxy_tmp,W,W_non_loc);

                    }
                }

                G2_step(h2G2xxxx_comm,psi_xxxx,pi_xxxx,K2xxxx,RK_result_G2xxxx,RK_weight);
                G2_step(h2G2xyxy_comm,psi_xyxy,pi_xyxy,K2xyxy,RK_result_G2xyxy,RK_weight);
                
            }
            rho_step(h1rho_comm,I,K1,RK_result_rho,RK_weight,Ns,Nb);

#else
            if(U!=0 && HF == false){

                thread th_h2G2xxxx_comm(commutator,ref(h2G2xxxx_comm),ref(G2xxxx_tmp),ref(h2));
                thread th_h2G2xyxy_comm(commutator,ref(h2G2xyxy_comm),ref(G2xyxy_tmp),ref(h2));
                thread th_coll_int(collision_int,ref(I),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W));

                th_h2G2xxxx_comm.join();
                th_h2G2xyxy_comm.join();
                th_coll_int.join();

                if(sparse==true){
                    if(SE=="scnd_brn"){
                        thread th_psi_sparse_xyxy(Psi_X_sparse,ref(psi_xyxy),ref(rho_tmp),ref(W));
                        thread th_psi_sparse_xxxx(Psi_X_sparse,ref(psi_xxxx),ref(rho_tmp),ref(W_non_loc));
                        th_psi_sparse_xyxy.join();
                        th_psi_sparse_xxxx.join();
                    }

                    if(SE=="GW"){
                        thread th_psi_sparse_xyxy(Psi_sparse,ref(psi_xyxy),ref(rho_tmp),ref(W));
                        thread th_psi_sparse_xxxx(Psi_sparse,ref(psi_xxxx),ref(rho_tmp),ref(W_non_loc));
                        thread th_pi_sparse_xxxx(Pixyxy_sparse,ref(pi_xxxx),ref(rho_tmp),ref(G2xyxy_tmp),ref(G2xxxx_tmp),ref(W),ref(W_non_loc));
                        thread th_pi_sparse_xyxy(Pixyxy_sparse,ref(pi_xyxy),ref(rho_tmp),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W),ref(W_non_loc));
                        th_pi_sparse_xxxx.join();
                        th_pi_sparse_xyxy.join();
                        th_psi_sparse_xyxy.join();
                        th_psi_sparse_xxxx.join();
                    }

                }
                else if(sparse==false){
                    if(SE=="scnd_brn"){
                        make_rho_les(rho_les,rho_tmp,Ns,Nb);
                        make_rho_grt(rho_grt,rho_tmp,Ns,Nb);

                        thread th_psi_xyxy(Psi_X,ref(psi_xyxy),ref(rho_les),ref(rho_grt),ref(W));
                        thread th_psi_xxxx(Psi_X,ref(psi_xxxx),ref(rho_les),ref(rho_grt),ref(W_non_loc));

                        th_psi_xxxx.join();
                        th_psi_xyxy.join();

                    }
                    else if(SE=="GW"){
                        make_rho_les(rho_les,rho_tmp,Ns,Nb);
                        make_rho_grt(rho_grt,rho_tmp,Ns,Nb);

                        thread th_psi_xyxy(Psi,ref(psi_xyxy),ref(rho_les),ref(rho_grt),ref(W));
                        thread th_psi_xxxx(Psi,ref(psi_xxxx),ref(rho_les),ref(rho_grt),ref(W_non_loc));
                        thread th_pi_xxxx(Pixyxy,ref(pi_xxxx),ref(rho_les),ref(rho_grt),ref(G2xyxy_tmp),ref(G2xxxx_tmp),ref(W),ref(W_non_loc));
                        thread th_pi_xyxy(Pixyxy,ref(pi_xyxy),ref(rho_les),ref(rho_grt),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W),ref(W_non_loc));

                        th_psi_xxxx.join();
                        th_psi_xyxy.join();
                        th_pi_xxxx.join();
                        th_pi_xyxy.join();
                    }
                }

                thread th_G2xxxx_step(G2_step,ref(h2G2xxxx_comm),ref(psi_xxxx),ref(pi_xxxx),ref(K2xxxx),ref(RK_result_G2xxxx),RK_weight);
                thread th_G2xyxy_step(G2_step,ref(h2G2xyxy_comm),ref(psi_xyxy),ref(pi_xyxy),ref(K2xyxy),ref(RK_result_G2xyxy),RK_weight);
                th_G2xxxx_step.join();
                th_G2xyxy_step.join();
            }
            rho_step(h1rho_comm,I,K1,RK_result_rho,RK_weight);

#endif
        }
        t+=dt;

        for(int i =0; i<Ns*Ns*Nb*Nb; ++i){
            rho[i] += RK_result_rho[i]*dt/6.0;
        }
        for(int i = 0; i<Ns2*Ns2*Nb2*Nb2; ++i){
            G2xxxx[i] += RK_result_G2xxxx[i]*dt/6.0;
            G2xyxy[i] += RK_result_G2xyxy[i]*dt/6.0;
        }
        
        cout<<double(t)/(tmax)*100<<"%\r";
        cout.flush();
//        progress.push_back(double(t)/(tmax)*100);
//        write_to_rfile("progress.txt",progress);
    }
    auto end_suc = chrono::steady_clock::now();
    cout<<"\n Execution Succesful, elapsed time:"<< chrono::duration<double>(end_suc - start).count()<< " sec \n";
    write_to_rfile("G2_max_eval.txt",evals_max);
    return 0;
}

//int RK4_DMD(vector<complex<double>> &rho,vector<vector<complex<double>>> &rho_diag,vector<vector<complex<double>>> &G2xyxy_diag,vector<vector<complex<double>>> &G2xxxx_diag,vector<complex<double>> &W,vector<complex<double>> &hsp,vector<vector<double>> &r,vector<vector<complex<double>>> &I_diag)
//{
//    cout<<"Beginning Time Stepping Procedure\n";
//
//    double tmax=double(time_steps)*dt_fixed;
//    auto start = chrono::steady_clock::now();
//    double RK_weight,dt,t,t_temp;
//    double eps=1e-3;
//
//    vector<complex<double>> I(Ns*Ns*Nb*Nb);
//    vector<complex<double>> h1rho_comm(Ns*Ns*Nb*Nb);
//    vector<complex<double>> rho_tmp(Ns*Ns*Nb*Nb);
//    vector<complex<double>> h(Ns*Ns*Nb*Nb);
//    vector<complex<double>> G2xxxx_tmp(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> G2xyxy_tmp(Ns2*Ns2*Nb2*Nb2);
//    cout<<G2xyxy_diag.size()<<"\n";
//    for(int t_index = 0; t_index < I_diag.size() - 1; ++t_index)
//    {
//        vector<double> progress;
////        if(t_index%5000 == 0 and t_index != 0){
////            if(HF == true){
////                write_to_cfile2D("rho_tot_HF.txt",rho_diag);
////            }
////            if(HF == false){
////                write_to_cfile2D("rho_tot_GKBA.txt",rho_diag);
////            }
////        }
////        for(int i = 0; i < Ns; ++i){
////            double num_dens = 0;
////            for(int a = 0; a < Nb; ++a){
////                num_dens+= rho[(i+a*Ns)*Nb*Ns + i + a*Ns].real();
////            }
////
////            if(num_dens>1.0 + eps || isnan(num_dens) || num_dens < -eps){
////                cout<<"ERROR:Unphysical Number Density on Site "<<(i+1)<<" Terminating Program. \nPercent Completed:"<<double(t)/(tmax)*100<<"% \n";
////                if(HF == true){
////                    write_to_cfile2D("rho_tot_HF.txt",rho_diag);
////                }
////                if(HF == false){
////                    write_to_cfile2D("rho_tot_GKBA.txt",rho_diag);
////                }
////
////                return 1;
////            }
////        }
//
//        vector<complex<double>> RK_result_rho(Ns*Ns*Nb*Nb);
//        vector<complex<double>> K1(Ns*Ns*Nb*Nb);
//
//
//        for(int k = 0; k<4; ++k){
//            vector<complex<double>> h1(Ns*Ns*Nb*Nb);
//
//            if(k == 1 or k == 2){
//                dt = dt_fixed/2;
//                RK_weight = 2;
//            }
//            else{
//                dt = dt_fixed;
//                RK_weight = 1;
//            }
//
//            if(k!=0){
//                t_temp = t + dt;
//            }
//            else{
//                t_temp = t;
//            }
//            for(int i = 0; i < rho.size(); ++i){
//                rho_tmp[i] = rho[i] + K1[i]*dt;
//            }
////            if(k == 0){
////                for(int i = 0; i < G2xxxx_tmp.size(); ++i){
////
//////                    G2xxxx_tmp[i] = 0.0*G2xyxy_diag[t_index][i];
////                    G2xyxy_tmp[i] = G2xyxy_diag[t_index][i];
////                }
////            }
////            else if(k == 1 or k == 2){
////                for(int i = 0; i < G2xxxx_tmp.size(); ++i){
////
//////                    G2xxxx_tmp[i] = 0.0*(G2xyxy_diag[t_index][i] + G2xyxy_diag[t_index+1][i])/2.0;
////                    G2xyxy_tmp[i] = (G2xyxy_diag[t_index+1][i] + G2xyxy_diag[t_index][i])/2.0;
////                }
////            }
////
////            else if(k == 3){
////                for(int i = 0; i < G2xxxx_tmp.size(); ++i){
////
//////                    G2xxxx_tmp[i] = 0.0*G2xyxy_diag[t_index+1][i];
////                    G2xyxy_tmp[i] = G2xyxy_diag[t_index+1][i];
////                }
////            }
//            if(k == 0){
//                for(int i = 0; i < I.size(); ++i){
//
////                    G2xxxx_tmp[i] = 0.0*G2xyxy_diag[t_index][i];
//                    I[i] = I_diag[t_index][i];
//                }
//            }
//            else if(k == 1 or k == 2){
//                for(int i = 0; i < I.size(); ++i){
//
////                    G2xxxx_tmp[i] = 0.0*(G2xyxy_diag[t_index][i] + G2xyxy_diag[t_index+1][i])/2.0;
//                    I[i] = (I_diag[t_index+1][i] + I_diag[t_index][i])/2.0;
//                }
//            }
//
//            else if(k == 3){
//                for(int i = 0; i < I.size(); ++i){
//
////                    G2xxxx_tmp[i] = 0.0*G2xyxy_diag[t_index+1][i];
//                    I[i] = I_diag[t_index+1][i];
//                }
//            }
//            h_HF(h1,hsp,rho_tmp,W,t_temp,Tp,wp,E,t0,r);
//
//
//            commutator(h1rho_comm,rho_tmp,h1);
////            collision_int(I,G2xxxx_tmp,G2xyxy_tmp,W);
//
//
//            rho_step(h1rho_comm,I,K1,RK_result_rho,RK_weight,Ns,Nb);
//
//        }
//
//        for(int i =0; i<Ns*Ns*Nb*Nb; ++i){
//            rho[i] += RK_result_rho[i]*dt/6.0;
//        }
//
//        rho_diag.push_back(rho);
//
//        t+=dt;
//        cout<<double(t)/(tmax)*100<<"%\r";
//        cout.flush();
//        progress.push_back(double(t)/(tmax)*100);
//    }
//
//
//
//    auto end_suc = chrono::steady_clock::now();
//    cout<<"\n Execution Succesful, elapsed time:"<< chrono::duration<double>(end_suc - start).count()<< " sec \n";
//    return 0;
//}




//int RK4_Hubbard_GKBA(vector<complex<double>> &rho,vector<vector<complex<double>>> &rho_diag)
//{
//    cout<<"Beginning Time Stepping Procedure\n";
//
//    double tmax=double(time_steps)*dt_fixed;
//    auto start = chrono::steady_clock::now();
//    double RK_weight,dt,t,t_temp;
//    dt=dt_fixed;
//    vector<complex<double>> I(Ns*Ns*Nb*Nb);
//
//
//    vector<complex<double>> W(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> U_t(Ns*Ns*Nb*Nb);
//    vector<complex<double>> U_t_dag(Ns*Ns*Nb*Nb);
//    vector<vector<complex<double>>> Y1;
//
//    vector<complex<double>> eye(Ns*Ns*Nb*Nb);
//    vector<vector<complex<double>>> I_tot_GKBA;
//    for(int i = 0; i <Ns; ++i){
//        eye[i*(Ns+1)]=1.0;
//    }
//    Y1.push_back(eye);
//
//    make_Hubbard_w(W,U);
//    for(int t_index = 0; t_index < time_steps; ++t_index)
//    {
//        rho_diag.push_back(rho);
//
//        for(int i = 0; i < Ns; ++i){
//            double num_dens = 0;
//            for(int a = 0; a < Nb; ++a){
//                num_dens+= rho[(i+a*Ns)*Nb*Ns + i + a*Ns].real();
//            }
//            if(num_dens > 1.05 || isnan(num_dens) || num_dens < -.05){
//                cout<<"ERROR:Unphysical Number Density on Site "<<(i+1)<<" Terminating Program. \nPercent Completed:"<<double(t)/(tmax)*100<<"% \n";
//                cout<<num_dens<<"\n";
//                write_to_cfile2D("I_GKBA.txt",I_tot_GKBA);
//
//                return 1;
//            }
//        }
//        vector<complex<double>> hHF(Ns*Ns*Nb*Nb);
//        vector<complex<double>> rho0(Ns*Ns*Nb*Nb);
//        for(int i = 0; i < rho.size(); ++i){
//            rho0[i] = rho[i];
//        }
//
//        for(int p = 0; p < 2; ++p){
//            vector<vector<complex<double>>> SE_lss_t2t1;
//            vector<vector<complex<double>>> SE_grt_t2t1;
//            vector<complex<double>> C(Ns*Ns*Nb*Nb);
//            vector<complex<double>> hHF_tmp(Ns*Ns*Nb*Nb);
//            vector<complex<double>> I_tmp(Ns*Ns*Nb*Nb);
//
//            vector<complex<double>> Cn(Ns*Ns*Nb*Nb);
//
//
//            h_HF_Hubbard(hHF_tmp,rho,W,t1,t2,t3,0,epsilon,0,t_temp, Tp, wp,E,t0);
//            scnd_brn_SE2(SE_lss_t2t1,SE_grt_t2t1,rho_diag,Y1,W);
//            collision_int_GKBA2(I_tmp,rho_diag,Y1,SE_lss_t2t1,SE_grt_t2t1,dt);
//            if(p==0){
//                for(int j = 0; j < hHF.size(); ++j){
//                    hHF[j]=hHF_tmp[j];
//                    I[j]=I_tmp[j];
//                }
//            }
//
//            else if(p==1){
//                for(int j = 0; j < hHF.size(); ++j){
//                    hHF[j]=(hHF[j]+hHF_tmp[j])/2.0;
//                    I[j]=(I[j]+I_tmp[j])/2.0;
//                }
//            }
//
//            U_delta(U_t,hHF);
//
//            for(int i = 0; i < Ns; ++i){
//                for(int j = 0;j <Ns;++j){
//                    U_t_dag[i*Ns + j] = conj(U_t[j*Ns + i]);
//                }
//            }
//
//            for(int i = 0; i < Ns; ++i){
//                I[i*(Ns+1)] = {0,I[i*(Ns+1)].imag()};
//            }
//
//            for(int j = 0; j < Ns; ++j){
//                for(int k = 0; k < Ns; ++k){
//                    Cn[j*Ns+k] = -im*dt*(I[j*Ns+k]+conj(I[k*Ns + j]));
//                }
//            }
//            for(int j = 0; j < C.size(); ++j){
//                C[j]+=Cn[j];
//            }
//
//
//            for(int i = 1; i < 3; ++i){
//                commutator(Cn,hHF,Cn);
//
//                for(int j = 0; j < C.size(); ++j){
//                    Cn[j]*=dt/(i+1);
//                }
//                for(int j = 0; j < C.size(); ++j){
//                    C[j]+=Cn[j];
//                }
//
//            }
//
//            for(int i = 0; i < rho.size(); ++i){
//                rho[i] = rho0[i] - im*C[i];
//            }
//
//            matrix_mult(rho,rho,U_t_dag);
//            matrix_mult(rho,U_t,rho);
//        }
//
//        h_HF_Hubbard(hHF,rho,W,t1,t2,t3,0,epsilon,0,t, Tp, wp,E,t0);
//
//        U_delta(U_t,hHF);
//
//        for(int i = 0; i < Ns; ++i){
//            for(int j = 0;j <Ns;++j){
//                U_t_dag[i*Ns + j] = conj(U_t[j*Ns + i]);
//            }
//        }
//        for(int i = 0; i < Y1.size(); ++i){
//            matrix_mult(Y1[i],Y1[i],U_t_dag);
//        }
//        Y1.push_back(eye);
//        I_tot_GKBA.push_back(I);
//
//        t+=dt;
//        cout<<double(t)/(tmax)*100<<"%\r";
//        cout.flush();
//
//
//    }
//
//    auto end_suc = chrono::steady_clock::now();
//    cout<<"\n Execution Succesful, elapsed time:"<< chrono::duration<double>(end_suc - start).count()<< " sec \n";
//    write_to_cfile2D("I_GKBA.txt",I_tot_GKBA);
//
//    return 0;
//}
//
//
//
//
//int RK4_GW(vector<complex<double>> &rho,vector<complex<double>> &rho_eq, vector<complex<double>> &G2xxxx, vector<complex<double>> &G2xyxy,vector<vector<complex<double>>> &rho_diag,vector<complex<double>> &W, vector<double> &GW_evals)
//{
//    cout<<"Beginning Time Stepping Procedure\n";
//
//    double tmax=double(time_steps)*dt_fixed;
//    auto start = chrono::steady_clock::now();
//    double RK_weight,dt,t,t_temp;
//    double eps=1e-3;
//    vector<complex<double>> I(Ns*Ns*Nb*Nb);
//
//    vector<complex<double>> rho_les(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> rho_grt(Ns2*Ns2*Nb2*Nb2);
//
//    vector<complex<double>> h2(Ns2*Ns2*Nb2*Nb2);
//
//    vector<complex<double>> W_non_loc(Ns2*Ns2*Nb2*Nb2);
//
//    vector<complex<double>> h2G2xxxx_comm(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> h2G2xyxy_comm(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> h1rho_comm(Ns*Ns*Nb*Nb);
//
//    vector<complex<double>> psi_xxxx(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> psi_xyxy(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> psi0_xxxx(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> psi0_xyxy(Ns2*Ns2*Nb2*Nb2);
//
//    vector<complex<double>> pi_xxxx(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> pi_xyxy(Ns2*Ns2*Nb2*Nb2);
//
//    vector<complex<double>> G2xxxx_tmp(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> G2xyxy_tmp(Ns2*Ns2*Nb2*Nb2);
//    vector<complex<double>> rho_tmp(Ns*Ns*Nb*Nb);
//
//    w_non_loc(W_non_loc,W);
//
//    if(sparse==true){
//        Psi_sparse(psi0_xyxy,rho_eq,W);
//        Psi_sparse(psi0_xxxx,rho_eq,W_non_loc);
//    }
//    else if(sparse==false){
//        make_rho_les(rho_les,rho_eq,Ns,Nb);
//        make_rho_grt(rho_grt,rho_eq,Ns,Nb);
//
//        Psi(psi0_xyxy,rho_les,rho_grt,W);
//        Psi(psi0_xxxx,rho_les,rho_grt,W_non_loc);
//    }
//
//
//    for(int t_index = 0; t_index < time_steps; ++t_index)
//    {
//        vector<double> progress;
//        if(t_index%5000 == 0 and t_index != 0){
//            if(HF == true){
//                write_to_cfile2D("rho_tot_HF.txt",rho_diag);
//            }
//            if(HF == false){
//                write_to_cfile2D("rho_tot_GKBA.txt",rho_diag);
//            }
//        }
//        for(int i = 0; i < Ns; ++i){
//            double num_dens = 0;
//            for(int a = 0; a < Nb; ++a){
//                num_dens+= rho[(i+a*Ns)*Nb*Ns + i + a*Ns].real();
//            }
//
//            if(num_dens>1.0 + eps || isnan(num_dens) || num_dens < -eps){
//                cout<<"ERROR:Unphysical Number Density on Site "<<(i+1)<<" Terminating Program. \nPercent Completed:"<<double(t)/(tmax)*100<<"% \n";
//                if(HF == true){
//                    write_to_cfile2D("rho_tot_HF.txt",rho_diag);
//                }
//                if(HF == false){
//                    write_to_cfile2D("rho_tot_GKBA.txt",rho_diag);
//                }
//
//                return 1;
//            }
//        }
//
//        vector<complex<double>> RK_result_rho(Ns*Ns*Nb*Nb);
//        vector<complex<double>> RK_result_G2xxxx(Ns2*Ns2*Nb2*Nb2);
//        vector<complex<double>> RK_result_G2xyxy(Ns2*Ns2*Nb2*Nb2);
//
//        vector<complex<double>> K1(Ns*Ns*Nb*Nb);
//        vector<complex<double>> K2xxxx(Ns2*Ns2*Nb2*Nb2);
//        vector<complex<double>> K2xyxy(Ns2*Ns2*Nb2*Nb2);
//
//        for(int k = 0; k<4; ++k){
//            vector<complex<double>> h1(Ns*Ns*Nb*Nb);
//
//            if(k == 1 or k == 2){
//                dt = dt_fixed/2;
//                RK_weight = 2;
//            }
//            else{
//                dt = dt_fixed;
//                RK_weight = 1;
//            }
//
//            if(k!=0){
//                t_temp = t + dt;
//            }
//            else{
//                t_temp = t;
//            }
//
//            for(int i = 0; i < Ns*Ns*Nb*Nb; ++i){
//                rho_tmp[i] = rho[i] + K1[i]*dt;
//            }
//
//            for(int i = 0; i < Ns2*Ns2*Nb2*Nb2; ++i){
//                G2xxxx_tmp[i] = G2xxxx[i] + K2xxxx[i]*dt;
//                G2xyxy_tmp[i] = G2xyxy[i] + K2xyxy[i]*dt;
//            }
//
//
//
//            h_HF_GW(h1,rho_tmp,rho_eq,WHF,GW_evals,t_temp,Tp,wp,E,t0);
//
//            HF2(h2,h1);
//            commutator(h1rho_comm,rho_tmp,h1);
//
//#if NOTHREADS
//
//            if(HF == false){
//
//                commutator(h2G2xxxx_comm,G2xxxx_tmp,h2);
//                commutator(h2G2xyxy_comm,G2xyxy_tmp,h2);
//                collision_int(I,G2xxxx_tmp,G2xyxy_tmp,W);
//
//                if(sparse==true){
//                    Psi_sparse(psi_xyxy,rho_tmp,W);
//                    Psi_sparse(psi_xxxx,rho_tmp,W_non_loc);
//                    Pixyxy_sparse(pi_xxxx,rho_tmp,G2xyxy_tmp,G2xxxx_tmp,W,W_non_loc);
//                    Pixyxy_sparse(pi_xyxy,rho_tmp,G2xxxx_tmp,G2xyxy_tmp,W,W_non_loc);
//                }
//                else if(sparse==false){
//                    make_rho_les(rho_les,rho_tmp,Ns,Nb);
//                    make_rho_grt(rho_grt,rho_tmp,Ns,Nb);
//
//                    Psi(psi_xyxy,rho_les,rho_grt,W);
//                    Psi(psi_xxxx,rho_les,rho_grt,W_non_loc);
//
//                    Pixyxy(pi_xxxx,rho_les,rho_grt,G2xyxy_tmp,G2xxxx_tmp,W,W_non_loc);
//                    Pixyxy(pi_xyxy,rho_les,rho_grt,G2xxxx_tmp,G2xyxy_tmp,W,W_non_loc);
//
//                }
//
//                G2_step_GW(h2G2xxxx_comm,psi_xxxx,psi0_xxxx,pi_xxxx,K2xxxx,RK_result_G2xxxx,RK_weight,Ns,Nb);
//                G2_step_GW(h2G2xyxy_comm,psi_xyxy,psi0_xyxy,pi_xyxy,K2xyxy,RK_result_G2xyxy,RK_weight,Ns,Nb);
//            }
//            rho_step(h1rho_comm,I,K1,RK_result_rho,RK_weight,Ns,Nb);
//#else
//            if(HF == false){
//
//                thread th_h2G2xxxx_comm(commutator,ref(h2G2xxxx_comm),ref(G2xxxx_tmp),ref(h2));
//                thread th_h2G2xyxy_comm(commutator,ref(h2G2xyxy_comm),ref(G2xyxy_tmp),ref(h2));
//                thread th_coll_int(collision_int,ref(I),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W));
//
//                th_h2G2xxxx_comm.join();
//                th_h2G2xyxy_comm.join();
//                th_coll_int.join();
//
//                if(sparse==true){
//                    thread th_psi_sparse_xyxy(Psi_sparse,ref(psi_xyxy),ref(rho_tmp),ref(W));
//                    thread th_psi_sparse_xxxx(Psi_sparse,ref(psi_xxxx),ref(rho_tmp),ref(W_non_loc));
//                    thread th_pi_sparse_xxxx(Pixyxy_sparse,ref(pi_xxxx),ref(rho_tmp),ref(G2xyxy_tmp),ref(G2xxxx_tmp),ref(W),ref(W_non_loc));
//                    thread th_pi_sparse_xyxy(Pixyxy_sparse,ref(pi_xyxy),ref(rho_tmp),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W),ref(W_non_loc));
//                    th_pi_sparse_xxxx.join();
//                    th_pi_sparse_xyxy.join();
//                    th_psi_sparse_xyxy.join();
//                    th_psi_sparse_xxxx.join();
//
//                }
//                else if(sparse==false){
//                    make_rho_les(rho_les,rho_tmp,Ns,Nb);
//                    make_rho_grt(rho_grt,rho_tmp,Ns,Nb);
//
//                    thread th_psi_xyxy(Psi,ref(psi_xyxy),ref(rho_les),ref(rho_grt),ref(W));
//                    thread th_psi_xxxx(Psi,ref(psi_xxxx),ref(rho_les),ref(rho_grt),ref(W_non_loc));
//                    thread th_pi_xxxx(Pixyxy,ref(pi_xxxx),ref(rho_les),ref(rho_grt),ref(G2xyxy_tmp),ref(G2xxxx_tmp),ref(W),ref(W_non_loc));
//                    thread th_pi_xyxy(Pixyxy,ref(pi_xyxy),ref(rho_les),ref(rho_grt),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W),ref(W_non_loc));
//
//                    th_psi_xxxx.join();
//                    th_psi_xyxy.join();
//                    th_pi_xxxx.join();
//                    th_pi_xyxy.join();
//                }
//
//
//                thread th_G2xxxx_step(G2_step_GW,ref(h2G2xxxx_comm),ref(psi_xxxx),ref(psi0_xxxx),ref(pi_xxxx),ref(K2xxxx),ref(RK_result_G2xxxx),RK_weight,Ns,Nb);
//                thread th_G2xyxy_step(G2_step_GW,ref(h2G2xyxy_comm),ref(psi_xyxy),ref(psi0_xyxy),ref(pi_xyxy),ref(K2xyxy),ref(RK_result_G2xyxy),RK_weight,Ns,Nb);
//                th_G2xxxx_step.join();
//                th_G2xyxy_step.join();
//            }
//            rho_step(h1rho_comm,I,K1,RK_result_rho,RK_weight,Ns,Nb);
//
//#endif
//        }
//
//        for(int i =0; i<Ns*Ns*Nb*Nb; ++i){
//            rho[i] += RK_result_rho[i]*dt/6.0;
//        }
//        for(int i = 0; i<Ns2*Ns2*Nb2*Nb2; ++i){
//            G2xxxx[i] += RK_result_G2xxxx[i]*dt/6.0;
//            G2xyxy[i] += RK_result_G2xyxy[i]*dt/6.0;
//        }
//        rho_diag.push_back(rho);
//
//        t+=dt;
//        cout<<double(t)/(tmax)*100<<"%\r";
//        cout.flush();
//        progress.push_back(double(t)/(tmax)*100);
//    }
//
//
//
//    auto end_suc = chrono::steady_clock::now();
//    cout<<"\n Execution Succesful, elapsed time:"<< chrono::duration<double>(end_suc - start).count()<< " sec \n";
//    return 0;
//}


#endif

 

