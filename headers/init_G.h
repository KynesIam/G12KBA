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

#ifndef init_G_h
#define init_G_h
#include "HF.h"
#include "clapack.h"
#include "linear_algebra.h"
//#include "mkl_lapacke.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>
using namespace std;

extern int Ns;
extern int Nb;
extern int Ns2;
extern int Nb2;



void rho0(vector<complex<double>> &rho,int p_number,vector<complex<double>> &hsp)
{
    vector<complex<double>> evecs(Ns*Nb*Ns*Nb);
    vector<double> evals(Ns*Nb);
    diagonalize(hsp,evecs,evals);
//    int degen=1;
//    vector<int> deg_store;
//    vector<int> idx_store;
//    int idx=0;
//    for(int i = 0; i < evals.size()-1; ++i){
//        if(evals[i] <= evals[i+1]+1e-10 and evals[i]>=evals[i+1]-1e-10){
//            degen+=1;
//
//        }
//        else{
//
//                deg_store.push_back(degen);
//                idx_store.push_back(idx);
//                idx+=degen;
//
//                degen=1;
//
//
//        }
//    }
//
    for(int i = 0; i<Ns*Nb; ++i){
        cout<<"Single particle energies"<<evals[i]<<"\n";
        for(int j = 0; j < Ns*Nb; ++j){
            rho[i*Ns*Nb + j] = 0;
            for(int k = 0;k<p_number;++k){
                if(k!=p_number-1){
                    rho[i*Ns*Nb + j] += (evecs[k*Ns*Nb + i]*conj(evecs[k*Ns*Nb+j]));
                }
                else if(k==p_number-1){
                    int count=0;
                    for(int l = k; l < evals.size(); ++l){
                        if(evals[k] <= evals[l]+1e-10 and evals[k] >= evals[l]-1e-10){
                            count+=1;
                        }
                        else{
                            break;
                        }
                    }
                    for(int l = k; l < k+count; ++l){
                        rho[i*Ns*Nb + j] += ((evecs[l*Ns*Nb + i]*conj(evecs[l*Ns*Nb+j])))/double(count);
                    }

                }
                
            }
           
        }
    }
}

//void rho0_tmp(vector<complex<double>> &rho,int p_number,vector<complex<double>> &hsp)
//{
//    vector<complex<double>> evecs(Ns*Nb*Ns*Nb);
//    vector<double> evals(Ns*Nb);
//    diagonalize(hsp,evecs,evals);
//    int degen=1;
//    vector<int> deg_store;
//    vector<int> idx_store;
//    int idx=0;
//    for(int i = 0; i < evals.size()-1; ++i){
//        if(evals[i] <= evals[i+1]+1e-10 and evals[i]>=evals[i+1]-1e-10){
//            degen+=1;
//
//        }
//        else{
//
////            if(degen>1){
//                deg_store.push_back(degen);
//                idx_store.push_back(idx);
//                idx+=degen;
//
//                degen=1;
////            }
//
//
//        }
//    }
//    int count2=0;
//    for(int i = 0; i < deg_store.size(); ++i){
//        count2+=deg_store[i];
//        if(count2 >= p_number){
//            for(int i = 0; i<Ns*Nb; ++i){
//                cout<<"Single particle energies"<<evals[i]<<"\n";
//                for(int j = 0; j < Ns*Nb; ++j){
//                    rho[i*Ns*Nb + j] = 0;
//                    int count=0;
//                    int k =0;
//                    while(count<p_number){
//                        for(int l = 0;l<idx_store.size();++l){
//                            if(k == idx_store[l]){
//                                for(int m = k; m < k+deg_store[l]; ++m){
//                                    rho[i*Ns*Nb + j] += (evecs[m*Ns*Nb + i]*conj(evecs[m*Ns*Nb+j]))/double(deg_store[l]);
//                                }
//                                k+=deg_store[l];
//
//                                break;
//
//                            }
//                        }
//                        count++;
//                    }
//
//                }
//            }
//            break;
//        }
//        else if(count2 == p_number){
//            vector<complex<double>> evecs(Ns*Nb*Ns*Nb);
//            vector<double> evals(Ns*Nb);
//            diagonalize(hsp,evecs,evals);
//            for(int i = 0; i<Ns*Nb; ++i){
//                cout<<"Single particle energies"<<evals[i]<<"\n";
//                for(int j = 0; j < Ns*Nb; ++j){
//                    rho[i*Ns*Nb + j] = 0;
//                    for(int k = 0;k<p_number;++k){
//                        rho[i*Ns*Nb + j] += (evecs[k*Ns*Nb + i]*conj(evecs[k*Ns*Nb+j]));
//                    }
//
//                }
//            }
//        }
//    }
//
//}

void rho0_GW(vector<complex<double>> &rho,int p_number)
{
    int N = int(sqrt(rho.size()));
    for(int i = 0; i < p_number; ++i){
        rho[i*(N+1)]=1.0;
    }
}

void init_G_test(vector<complex<double>> &rho)
{
    for(int i = 0; i < int(Ns/2); ++i){
        for(int a = 0; a < int(Nb); ++a){
            rho[(i+a*Ns)*Nb*Ns + i + a*Ns] = 1;
        }
    }
}




#endif





