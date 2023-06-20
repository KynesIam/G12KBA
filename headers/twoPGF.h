//
//  twoPGF.h
//  
//
//  Created by Cian Reeves on 6/5/23.
//

#ifndef twoPGF_h
#define twoPGF_h
#include "print_matrix.h"
#include <vector>
#include <iostream>
#include <complex>
#include <string>
#include <fstream>
using namespace std;


void two_PGF(vector<complex<double>> &G2_full,vector<complex<double>> &G2_corr,vector<complex<double>> &G1){
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
                                    int idx1=IDX_8D(ks,ls,kb,lb,is,js,ib,jb);

                                    int idx2=IDX_4D(ls,lb,js,jb);
                                    int idx3=IDX_4D(ks,kb,is,ib);
                                    G2_full[idx1] = G2_corr[idx1] - G1[idx3]*G1[idx2];
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
}
#endif /* twoPGF_h */
