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


#ifndef collision_integral_h
#define collision_integral_h
#include "kron_delta.h"
#include "print_matrix.h"
#include <vector>
#include <complex>
using namespace std;
extern int Ns;
extern int Nb;
extern double dt_fixed;
void collision_int(vector<complex<double>> &I, vector<complex<double>> &G2xxxx,vector<complex<double>> &G2xyxy,vector<complex<double>> &W)
{
    /*
     Function to calculate the collision integral in the G1-G2 formulation for a given interaction matrix and 2PGF
     Args:
        I: Ns*Ns*Nb*Nb vector to store the the collision integral
        G2xxxx: (Ns*Ns*Nb*Nb)**2 vector that holds the current value of the 2PGF with both particles in same spin sector
        G2xyxy: (Ns*Ns*Nb*Nb)**2 vector that holds the current value of the 2PGF with both particles in opposite spin sector
        W: (Ns*Ns*Nb*Nb)**2 vector that holds the two-body interaction term
        
     */
    int Ns2 = int(pow(Ns,2));
    int Nb2 = int(pow(Nb,2));
    complex<double> im = {0,1};
    fill(I.begin(),I.end(),0);
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
#if NOTHREADS
    #pragma omp parallel for collapse(10)
#endif
    for(int ls = 0; ls < Ns; ++ls){
        for(int lb = 0; lb < Nb; ++lb){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    for(int is = 0; is < Ns; ++is){
                        for(int ib = 0; ib < Nb; ++ib){
                            for(int ms = 0; ms < Ns; ++ms){
                                for(int mb = 0; mb < Nb; ++mb){
                                    for(int ns = 0; ns < Ns; ++ns){
                                        for(int nb = 0; nb < Nb; ++nb){
                                            int idx1 = IDX_4D(js,jb,ls,lb);
                                            int idx2 = IDX_8D(is,js,ib,jb,ns,ms,nb,mb);
                                            int idx3 = IDX_8D(ns,ms,nb,mb,is,ls,ib,lb);
                                            if(is==ls && is==ns && is==ms && ib==lb && ib==nb && ib==mb){
                                                I[idx1] += -im*G2xyxy[idx2]*W[idx3];
                                            }
                                            else{
                                                I[idx1] += -im*(G2xxxx[idx2]+G2xyxy[idx2])*W[idx3];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


#endif

