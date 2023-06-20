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

#ifndef HF_h
#define HF_h
#include "kron_delta.h"
#include "print_matrix.h"
#include <vector>
#include <iostream>
#include <complex>
#include <string>
#include <fstream>
using namespace std;
extern int Ns;
extern int Nb;
extern int Ns2;
extern int Nb2;
extern double a;
extern double offset;
extern double alpha;
extern bool cluster;
extern double lambda;
extern bool pb;
extern vector<double> epsilon;
#define PI 3.14159265359

double Heaviside(double x)
{
    if(x<0.0){
        return 0.0;
    }
    if(x>=0.0){
        return 1.0;
    }
    return 0;
}

double distance(vector<double> &r1, vector<double> &r2){
    double d = 0;
    for(int i = 0; i < r1.size(); ++i){
        d += pow(abs(r1[i]-r2[i]),2);
    }
    return sqrt(d);
}

double pulse(double t, double Tp, double wp,double E,double t0)
{
    return E*exp(-(pow((t-t0)/Tp,2))/2.0);
//    return E*sin(t);
//    return E*Heaviside(1-abs(1-2*(t-t0)/Tp))*sin(wp*(t-t0))*pow(sin(PI*(t-t0)/Tp),2);
}




void h_sp_Hubbard_nd(vector<complex<double> > &h,vector<vector<double>> &r,  double t1)
{
    /*
     Function to initiate the time independent part of the single particle Hamiltonian for a lattice model with lattice points defined by vectors r_1,...,r_N
     Args:
        h: Ns*Ns*Nb*Nb vector that stores the single particle Hamiltonian
        r: vector of vectors to define the lattice points
        t1: Hopping strength between sites(Hoppings between orbitals can easily be implemented too)
     */

#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
    double d;
    double dmax=distance(r[0],r[r.size()-1]);
    if(pb==true){
        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        int idx1 = IDX_4D(js,jb,is,ib);
                        d=distance(r[is],r[js]);
                        
                        if(js != is and ib == jb and cluster == true){
                            d = min(d,dmax+a-d);
                            h[idx1] += -t1*exp(-alpha*(d-a));
                        }
                        
                        if(js == is and ib == jb){
                            h[idx1]+=epsilon[jb];
                        }
                        else if(cluster == false and ((d<= 1.0 + 1e-10 and d >= 1.0 - 1e-10) or (d <= Ns-1 + 1e-10 and d >= Ns-1 - 1e-10)) and ib==jb and ib == 0){
                            h[idx1] += -t1;
                        }
//                        else if(cluster == false and ((d<= 1.0 + 1e-10 and d >= 1.0 - 1e-10) or (d <= Ns-1 + 1e-10 and d >= Ns-1 - 1e-10)) and ib==jb and ib == 1){
//                            h[idx1] += t1;
//                        }
                    }
                }
            }
        }
    }
    
    else if(pb==false){
        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        int idx1 = IDX_4D(js,jb,is,ib);
                        d = distance(r[is],r[js]);
                        
                        if(js != is and cluster == true){
                            h[idx1] = -t1*exp(-alpha*(d-a));
                        }
                        
                        else if(cluster == false and  (d <= 1.0 +1e-10 and d >= 1.0 - 1e-10)){
                            h[idx1] = -t1;
                        }
                    }
                }
            }
        }
    }
}






void h_HF(vector<complex<double>> &h_HF,vector<complex<double>> &h_sp,vector<complex<double>> &rho, vector<complex<double>> &W, double t, double Tp, double wp,double E,double t0,vector<vector<double>> &r)
{
    /*
     Function that takes the single particle Hamiltonian and includes the Hartree-Fock contribution as well as the time dependent perturbation
     Args:
        h_HF: Ns*Ns*Nb*Nb vector that stores the full Hartree-Fock Hamiltonian
        h_sp: Ns*Ns*Nb*Nb vector that stores the single particle Hamiltonian
        rho: Ns*Ns*Nb*Nb vector that stores the density matrix at a given time
        W: (Ns*Ns*Ns*Nb)**2 vector that stores the full two-body interaction matrix
        t: time at which the HF Hamiltonian is being computed
        Tp: Parameter that determines temporal width of time dependent field(if applicable)
        wp: Parameter that determines osscilation frequency of time dependent field(if applicable)
        E: Paramter that determines strength of time dependent field
        t0: Parameter that determines midpoint or start of time dependent field(depending on perturbation type
        r: vector of vectors to define the lattice points
        
     */
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    
                    int idx1 = IDX_4D(js,jb,is,ib);
                    h_HF[idx1] = h_sp[idx1];

                    if(js==is){

                        if(dim == 3){

                            h_HF[idx1] += (r[is][0] + r[is][1] + r[is][2])*pulse(t, Tp, wp,E,t0)/sqrt(3);
                        }
                        else if(dim == 1){
//                        h_HF[idx1] += -E*cos(PI*(r[is][0] + t - t0)/2)*Heaviside(t-(t0+r[is][0]));
                            h_HF[idx1] += r[is][0]*pulse(t, Tp, wp,E,t0);
//                        h_HF[idx1]+= -cos(PI*r[is][0]/2.0)*pulse(t, Tp, wp,E,t0);
//                        h_HF[idx1] += -cos(2.0*PI*(r[is][0] - (t-t0))/lambda)*pulse(r[is][0] - (t-t0), Tp, wp,E,0);
//                        h_HF[idx1]+= -cos(2.0*PI*r[is][0]/4.0)*pulse(t, Tp, wp,E,t0);
                        }

                    }
                    for(int ks = 0; ks < Ns; ++ks){
                        for(int kb = 0; kb < Nb; ++kb){
                            for(int ls = 0; ls < Ns; ++ls){
                                for(int lb = 0; lb < Nb; ++lb){
                                    int idx2 = IDX_4D(ks,kb,ls,lb);
                                    int idx3 = IDX_8D(js,is,jb,ib,ls,ks,lb,kb);
                                    int idx4 = IDX_8D(ls,is,lb,ib,js,ks,jb,kb);
                                    if(is==js && is==ks && is==ls && ib==lb && ib==kb && ib==lb){
                                        h_HF[idx1] += (W[idx3])*(rho[idx2]);
                                    }
                                    else{
                                        h_HF[idx1] += (2.0*W[idx3] - W[idx4])*(rho[idx2]);
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

void HF2(vector<complex<double>> &h2, vector<complex<double>> &h)
{
  
    /*
     Function to compute two particle Hartree-Fock Hamiltonian given the single particle HF Hamiltonian
     Args:
        h2: (Ns*Ns*Nb*Nb)**2 vector to store two particle HF Hamiltonian
        h: (Ns*Ns*Nb*Nb) vector that stores single particle HF Hamiltonian
     */
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
                                    
                                    h2[idx1] = h[idx2]*kron_delta(ns,is)*kron_delta(nb,ib) - h[idx3]*kron_delta(ms,ls)*kron_delta(mb,lb);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



//void h_HF_GW(vector<complex<double>> &h_HF,vector<complex<double>> &rho,vector<complex<double>> &rho_eq, vector<complex<double>> &W, vector<double> &GW_evals,double t, double Tp, double wp,double E,double t0)
//{
//#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
//#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
//    for(int is = 0; is < Ns; ++is){
//        for(int ib = 0; ib < Nb; ++ib){
//            for(int js = 0; js < Ns; ++js){
//                for(int jb = 0; jb < Nb; ++jb){
//
//                    int idx1 = IDX_4D(js,jb,is,ib);
//                    if(is == js){
//                        h_HF[idx1] = GW_evals[is];
//                    }
//
//                    if(is != js){
////                        h_HF[idx1] += pulse(t,Tp,wp,E,t0);
//                    }
//
//                    for(int ks = 0; ks < Ns; ++ks){
//                        for(int kb = 0; kb < Nb; ++kb){
//                            for(int ls = 0; ls < Ns; ++ls){
//                                for(int lb = 0; lb < Nb; ++lb){
//                                    int idx2 = IDX_4D(ks,kb,ls,lb);
//                                    int idx3 = IDX_8D(js,is,jb,ib,ls,ks,lb,kb);
//                                    int idx4 = IDX_8D(ls,is,lb,ib,js,ks,jb,kb);
//                                    if(is==js && is==ks && is==ls && ib==lb && ib==kb && ib==lb){
//                                        h_HF[idx1] += (W[idx3])*(rho[idx2]-rho_eq[idx2]);
//                                    }
//                                    else{
//                                        h_HF[idx1] += (2.0*W[idx3] - W[idx4])*(rho[idx2]-rho_eq[idx2]);
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//}
#endif
