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

#ifndef two_body_h
#define two_body_h

#include "kron_delta.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
extern double decay_rate;
extern double coulomb_scaling;
extern string decay_type;
//void make_Hubbard_w(vector<complex<double>> &W,double U)
//{
//#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
//
//    for(int is = 0; is < Ns; ++is){
//        for(int ib = 0; ib < Nb; ++ib){
//            for(int ls = 0; ls < Ns; ++ls){
//                for(int lb = 0; lb < Nb; ++lb){
//                    for(int ms = 0; ms < Ns; ++ms){
//                        for(int mb = 0; mb < Nb; ++mb){
//                            for(int ns = 0; ns < Ns; ++ns){
//                                for(int nb = 0; nb < Nb; ++nb){
//                                    
//                                    int idx1 = IDX_8D(ns,is,nb,ib,ms,ls,mb,lb);
//                                    
//                                    W[idx1] =  U*kron_delta(is,ls)*kron_delta(is,ms)*kron_delta(is,ns)*kron_delta(ib,nb)*kron_delta(lb,mb)*(kron_delta(ib,lb));
//                                    W[idx1] += U*kron_delta(is,ls)*kron_delta(is,ms)*kron_delta(is,ns)*kron_delta(ib,nb)*kron_delta(lb,mb)*(1.0-kron_delta(ib,lb));
//
//                                    if(EHM==true){
//                                        for(int n = 1; n < Ns; ++n){
//                                            if(abs(ls-ns)==n){
//                                                if(decay_type == "exp"){
//                                                    W[idx1] += U*exp(-decay_rate*n)*kron_delta(is,ns)*kron_delta(ls,ms)*kron_delta(ib,nb)*kron_delta(lb,mb);
//                                                }
//                                                else if(decay_type == "coulomb"){
//                                                    W[idx1] += (U/double(n))*coulomb_scaling*kron_delta(is,ns)*kron_delta(ls,ms)*kron_delta(ib,nb)*kron_delta(lb,mb);
//
//                                                }
//                                            }
//                                        }
//                                    }
////                                    if(is!=ns and ls!=ms){
////                                        W[idx1]+=.001*U;
////                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//}
void make_Hubbard_w(vector<complex<double>> &W,double U,vector<vector<double>> r)
{
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)

    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int ls = 0; ls < Ns; ++ls){
                for(int lb = 0; lb < Nb; ++lb){
                    for(int ms = 0; ms < Ns; ++ms){
                        for(int mb = 0; mb < Nb; ++mb){
                            for(int ns = 0; ns < Ns; ++ns){
                                for(int nb = 0; nb < Nb; ++nb){
                                    
                                    int idx1 = IDX_8D(ns,is,nb,ib,ms,ls,mb,lb);
                                    
                                    W[idx1] =  U*kron_delta(is,ls)*kron_delta(is,ms)*kron_delta(is,ns)*kron_delta(ib,nb)*kron_delta(lb,mb)*(kron_delta(ib,lb));
                                    W[idx1] += U*kron_delta(is,ls)*kron_delta(is,ms)*kron_delta(is,ns)*kron_delta(ib,nb)*kron_delta(lb,mb)*(1.0-kron_delta(ib,lb));

                                    if(EHM==true){
                                        if(is != ls){
                                            double d = distance(r[is],r[ls]);
                                            if(decay_type == "exp"){
                                                W[idx1] += U*exp(-decay_rate*d)*kron_delta(is,ns)*kron_delta(ls,ms)*kron_delta(ib,nb)*kron_delta(lb,mb);
                                            }
                                            else if(decay_type == "coulomb"){
                                                W[idx1] += (U/d)*coulomb_scaling*kron_delta(is,ns)*kron_delta(ls,ms)*kron_delta(ib,nb)*kron_delta(lb,mb);
                                                
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


void w_non_loc(vector<complex<double>> &W_non_loc,vector<complex<double>> &W)
{
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int ls = 0; ls < Ns; ++ls){
                for(int lb = 0; lb < Nb; ++lb){
                    for(int ms = 0; ms < Ns; ++ms){
                        for(int mb = 0; mb < Nb; ++mb){
                            for(int ns = 0; ns < Ns; ++ns){
                                for(int nb = 0; nb < Nb; ++nb){
                                    if(is==ls && is==ms && is==ns && ib==lb && ib==mb && ib==nb){
                                        int idx1 = IDX_8D(ns,is,nb,ib,ms,ls,mb,lb);
                                        W_non_loc[idx1] = 0;
                                    }
                                    else{
                                        int idx1 = IDX_8D(ns,is,nb,ib,ms,ls,mb,lb);
                                        W_non_loc[idx1] = W[idx1];
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

void scaled_W(vector<complex<double>> &scaled_W,vector<complex<double>> &W,double scaling_factor)
{
    for(int i = 0; i < W.size(); ++i){
        scaled_W[i] = scaling_factor*W[i];
    }
}





#endif
