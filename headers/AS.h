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

#ifndef AS_h
#define AS_h


#include "HF.h" //Holds Hartree-Fock Hamiltonian
#include "psi.h"//Holds two particle occupation term
#include "pi.h" //Function for the polarization term in GW
#include "G1_G2_t_step.h"//Holds update steps for the one and two particle Green's function
#include "linear_algebra.h"//Commutators for one and two particle Green's function with respective HF hamiltonian
#include "collision_integral.h"//Collision integral, describes many body effects on one particle Green's function
#include "fermi_func.h"//Fermi function for AS procedure
#include "read_input.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <thread>
#include <complex>
#include <fstream>
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
extern double dt_fixed;
extern double quench_strength;
extern int Nq;
extern int Ns;
extern int Ns2;
extern int Nb;
extern int Nb2;
extern string q_type;
extern vector<double> epsilon;
extern double AS_rate;
extern double AS_midpoint;
extern vector<vector<double>> r;

int Adiabatic_switching(vector<complex<double> > &rho, vector<complex<double>> &G2xxxx,vector<complex<double> > &G2xyxy,vector<complex<double>> &W,vector<complex<double>> &hsp, vector<vector<complex<double>>> &rho_diag)
{
    /*
     Function that performs adiabatic switching procedure to prepare GKBA ground state.  Achieved by time evolving as the entire W matrix is scaled by the fermi-dirac distribution as a function
     of the evolution time
     Args:
        I: Ns*Ns*Nb*Nb vector that stores the non-interacting GS density matrix and at the end stores the GKBA ground state density matrix
        G2xxxx: (Ns*Ns*Nb*Nb)**2 vector that holds the current value of the 2PGF with both particles in same spin sector
        G2xyxy: (Ns*Ns*Nb*Nb)**2 vector that holds the current value of the 2PGF with both particles in opposite spin sector
        W: (Ns*Ns*Nb*Nb)**2 vector that holds the two-body interaction term
        hsp: Ns*Ns*Nb*Nb vector that stores time independent part of single particle Hamiltonian
        rho_diag: n_steps*(Ns*Ns*Nb*Nb) vector that stores the density matrix over the entire AS time period
     */
    cout<<"Beginning Adiabatic Switching Procedure\n";
    auto start = chrono::steady_clock::now();
    double RK_weight,t_RK,t,dt;
    
    double scaling_factor;
    
    vector<complex<double>> I(Ns*Ns*Nb*Nb);
    
    vector<complex<double>> rho_les(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> rho_grt(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> h2(Ns2*Ns2*Nb2*Nb2);
    
    vector<complex<double>> W_AS(Ns2*Ns2*Nb2*Nb2);
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
    
    double eps=1e-3;
    while(abs(scaling_factor-1)<=1e-10)
    {
        rho_diag.push_back(rho);

        for(int i = 0; i < Ns; ++i){
            //Checks if any particle densities on a given site are unphysical, either negative or greater than the number of orbitals per site
            double num_dens = 0;
            for(int a = 0; a < Nb; ++a){
                num_dens += rho[(i+a*Ns)*Nb*Ns + i + a*Ns].real();
            }
            
            if(num_dens> Nb+ eps || isnan(num_dens) || num_dens < -eps){
                cout<<"ERROR:Unphysical Number Density on Site "<<(i+1)<<"\n";
                if(HF == true){
                    write_to_cfile2D("rho_tot_HF.txt",rho_diag);
                }
                if(HF == false){
                    write_to_cfile2D("rho_tot_GKBA.txt",rho_diag);
                }

                return 1;
            }
        }
        vector<complex<double>> RK_result_rho(Ns*Ns*Nb*Nb);
        vector<complex<double>> K1(Ns*Ns*Nb*Nb);

        vector<complex<double>> RK_result_G2xxxx(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> RK_result_G2xyxy(Ns2*Ns2*Nb2*Nb2);

        vector<complex<double>> K2xxxx(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> K2xyxy(Ns2*Ns2*Nb2*Nb2);
        t_RK = t;
        
        for(int k = 0; k<4; ++k)
        {
            
            vector<complex<double>> h1(Ns*Ns*Nb*Nb);
            if(k == 1 or k == 2){
                dt = .01/2;
                RK_weight = 2;
            }
            else{
                dt = .01;
                RK_weight = 1;
            }
            if(k!=0){
                t_RK = t+dt;
            }
            else{
                t_RK = t;
            }

            for(int i = 0; i < Ns*Ns*Nb*Nb; ++i){
                rho_tmp[i] = rho[i] + K1[i]*dt;
            }

            for(int i = 0; i < Ns2*Ns2*Nb2*Nb2; ++i){
                G2xxxx_tmp[i] = G2xxxx[i] + K2xxxx[i]*dt;
                G2xyxy_tmp[i] = G2xyxy[i] + K2xyxy[i]*dt;
            }
            
            scaling_factor=fermi(t_RK,AS_midpoint,AS_rate);
            scaled_W(W_AS,W,scaling_factor);
            
            w_non_loc(W_non_loc,W_AS);
            h_HF(h1,hsp,rho_tmp,W_AS,1.0,1.0,1.0,0.0,1.0,r);
            HF2(h2,h1);
            commutator(h1rho_comm,rho_tmp,h1);
#if NOTHREADS
            if(U!=0 && HF == false){

                commutator(h2G2xxxx_comm,G2xxxx_tmp,h2);
                commutator(h2G2xyxy_comm,G2xyxy_tmp,h2);
                collision_int(I,G2xxxx_tmp,G2xyxy_tmp,W_AS);


                if(sparse==true){
                    if(SE=="scnd_brn"){
                        Psi_X_sparse(psi_xyxy,rho_tmp,W_AS);
                        Psi_X_sparse(psi_xxxx,rho_tmp,W_non_loc);
                        
                    }
                    
                    if(SE=="GW"){
                        Psi_sparse(psi_xyxy,rho_tmp,W_AS);
                        Psi_sparse(psi_xxxx,rho_tmp,W_non_loc);
                        Pixyxy_sparse(pi_xxxx,rho_tmp,G2xyxy_tmp,G2xxxx_tmp,W_AS,W_non_loc);
                        Pixyxy_sparse(pi_xyxy,rho_tmp,G2xxxx_tmp,G2xyxy_tmp,W_AS,W_non_loc);
                    }
                    
                }
                else if(sparse==false){
                    if(SE=="scnd_brn"){
                        make_rho_les(rho_les,rho_tmp,Ns,Nb);
                        make_rho_grt(rho_grt,rho_tmp,Ns,Nb);
                        
                        Psi_X(psi_xyxy,rho_les,rho_grt,W_AS);
                        Psi_X(psi_xxxx,rho_les,rho_grt,W_non_loc);

                    }
                    else if(SE=="GW"){
                        make_rho_les(rho_les,rho_tmp,Ns,Nb);
                        make_rho_grt(rho_grt,rho_tmp,Ns,Nb);
                        
                        Psi(psi_xyxy,rho_les,rho_grt,W_AS);
                        Psi(psi_xxxxrho_les,rho_grt,W_non_loc);
                        Pixyxy(pi_xxxx,rho_les,rho_grt,G2xyxy_tmp,G2xxxx_tmp,W_AS,W_non_loc);
                        Pixyxy(pi_xyxy,rho_les,rho_grt,G2xxxx_tmp,G2xyxy_tmp,W_AS,W_non_loc);

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
                thread th_coll_int(collision_int,ref(I),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W_AS));

                th_h2G2xxxx_comm.join();
                th_h2G2xyxy_comm.join();
                th_coll_int.join();

                if(sparse==true){
                    if(SE=="scnd_brn"){
                        thread th_psi_sparse_xyxy(Psi_X_sparse,ref(psi_xyxy),ref(rho_tmp),ref(W_AS));
                        thread th_psi_sparse_xxxx(Psi_X_sparse,ref(psi_xxxx),ref(rho_tmp),ref(W_non_loc));
                        th_psi_sparse_xyxy.join();
                        th_psi_sparse_xxxx.join();
                    }
                    
                    if(SE=="GW"){
                        thread th_psi_sparse_xyxy(Psi_sparse,ref(psi_xyxy),ref(rho_tmp),ref(W_AS));
                        thread th_psi_sparse_xxxx(Psi_sparse,ref(psi_xxxx),ref(rho_tmp),ref(W_non_loc));
                        thread th_pi_sparse_xxxx(Pixyxy_sparse,ref(pi_xxxx),ref(rho_tmp),ref(G2xyxy_tmp),ref(G2xxxx_tmp),ref(W_AS),ref(W_non_loc));
                        thread th_pi_sparse_xyxy(Pixyxy_sparse,ref(pi_xyxy),ref(rho_tmp),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W_AS),ref(W_non_loc));
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
                        
                        thread th_psi_xyxy(Psi_X,ref(psi_xyxy),ref(rho_les),ref(rho_grt),ref(W_AS));
                        thread th_psi_xxxx(Psi_X,ref(psi_xxxx),ref(rho_les),ref(rho_grt),ref(W_non_loc));

                        th_psi_xxxx.join();
                        th_psi_xyxy.join();

                    }
                    else if(SE=="GW"){
                        make_rho_les(rho_les,rho_tmp,Ns,Nb);
                        make_rho_grt(rho_grt,rho_tmp,Ns,Nb);
                        
                        thread th_psi_xyxy(Psi,ref(psi_xyxy),ref(rho_les),ref(rho_grt),ref(W_AS));
                        thread th_psi_xxxx(Psi,ref(psi_xxxx),ref(rho_les),ref(rho_grt),ref(W_non_loc));
                        thread th_pi_xxxx(Pixyxy,ref(pi_xxxx),ref(rho_les),ref(rho_grt),ref(G2xyxy_tmp),ref(G2xxxx_tmp),ref(W_AS),ref(W_non_loc));
                        thread th_pi_xyxy(Pixyxy,ref(pi_xyxy),ref(rho_les),ref(rho_grt),ref(G2xxxx_tmp),ref(G2xyxy_tmp),ref(W_AS),ref(W_non_loc));

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

        for(int i =0; i<Ns*Ns*Nb*Nb;++i){
            rho[i] += RK_result_rho[i]*dt/6.0;
        }
        for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i){
            G2xxxx[i] += RK_result_G2xxxx[i]*dt/6.0;
            G2xyxy[i] += RK_result_G2xyxy[i]*dt/6.0;
        }
        t+=dt;
        cout<<"Adiabatic switching"<<scaling_factor*100<<"%"<<"complete"<<"\r";
        cout.flush();

    }
    auto end = chrono::steady_clock::now();
    cout<<"Adiabatic switching complete, elapsed time:"<< chrono::duration<double>(end - start).count()<< " sec \n";
    return 0;
}


#endif 
