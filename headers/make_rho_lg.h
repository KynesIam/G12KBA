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


#ifndef make_rho_lg_h
#define make_rho_lg_h
#include <vector>
#include <complex>

void make_rho_grt(vector<complex<double>> &rho_grt, vector<complex<double>> &rho, int Ns, int Nb)
{
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
#if NOTHREADS
    #pragma omp parallel for collapse(8)
#endif
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int ls = 0; ls < Ns; ++ls){
                for(int lb = 0; lb < Nb; ++lb){
                    for(int ms = 0; ms < Ns; ++ms){
                        for(int mb = 0; mb < Nb; ++mb){
                            for(int ns = 0; ns < Ns; ++ns){
                                for(int nb = 0; nb < Nb; ++nb){
                                    
                                    int idx1 = IDX_8D(ns,ms,nb,mb,is,ls,ib,lb);
                                    int idx2 = IDX_4D(ms,mb,ls,lb);
                                    int idx3 = IDX_4D(is,ib,ns,nb);
                                    
                                    rho_grt[idx1] = rho[idx3]*(rho[idx2] - kron_delta(ls,ms)*kron_delta(lb,mb));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void make_rho_les(vector<complex<double>> &rho_les, vector<complex<double>> &rho, int Ns, int Nb)
{
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
#if NOTHREADS
    #pragma omp parallel for collapse(8)
#endif
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int ls = 0; ls < Ns; ++ls){
                for(int lb = 0; lb < Nb; ++lb){
                    for(int ms = 0; ms < Ns; ++ms){
                        for(int mb = 0; mb < Nb; ++mb){
                            for(int ns = 0; ns < Ns; ++ns){
                                for(int nb = 0; nb < Nb; ++nb){
                                    
                                    int idx1 = IDX_8D(ns,ms,nb,mb,is,ls,ib,lb);
                                    int idx2 = IDX_4D(ms,mb,ls,lb);
                                    int idx3 = IDX_4D(is,ib,ns,nb);
                                    
                                    rho_les[idx1] = rho[idx2]*(rho[idx3] - kron_delta(ns,is)*kron_delta(nb,ib));
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
