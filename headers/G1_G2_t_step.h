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


#ifndef G1_G2_tstep_h
#define G1_G2_tstep_h
#include <vector>
#include <complex>
using namespace std;
extern int Ns;
extern int Nb;
void rho_step(vector<complex<double>> &h1rho_comm, vector<complex<double>> &I,vector<complex<double>> &K,vector<complex<double>> &RK_result, double w)
{
    /*
     Function that performs one update step in the Runge-Kutta time stepping procedure and stores the neccessary data for the full algorithm
     Args:
        h1rho_comm: Ns*Ns*Nb*Nb vector that stores the commutator of the Hartree-Fock Hamiltonian and the density matrix
        I: Ns*Ns*Nb*Nb vector that stores the collision integral computed at each time step
        K: Ns*Ns*Nb*Nb vector that stores the current update step in the RK procedure to be used at the next step of the algorithm
        RK_result: Ns*Ns*Nb*Nb vector that stores the appropriately weighted sum of all previous RK steps for use in the final update step
        w: determines the weight that each K contributes to the final result
     
     */
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
    complex<double> im = {0,1};
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                   int idx1 = IDX_4D(js,jb,is,ib);
                   int idx2 = IDX_4D(is,ib,js,jb);
    
                   K[idx1] = -im*h1rho_comm[idx1] - (I[idx1] + conj(I[idx2]));
                   RK_result[idx1] += w*K[idx1];
               }
           }
       }
   }
}


void G2_step(vector<complex<double>> &h2G2_comm, vector<complex<double>> &psi, vector<complex<double>> &pi, vector<complex<double>> &K,vector<complex<double>> &RK_result, double w)
{
    /*
     Function that performs one update step in the Runge-Kutta time stepping procedure and stores the neccessary data for the full algorithm
     Args:
        h2G2_comm: (Ns*Ns*Nb*Nb)**2 vector that stores the commutator of the 2 particle Hartree-Fock Hamiltonian and the 2 particle Green's function
        psi: (Ns*Ns*Nb*Nb)**2 vector that stores the contribution to correlations from two particle scattering events
        pi: (Ns*Ns*Nb*Nb)**2 vector that stores the contribution to correlations from the polarization effects present in the GW approximation
        K: Ns*Ns*Nb*Nb vector that stores the current update step in the RK procedure to be used at the next step of the algorithm
        RK_result: Ns*Ns*Nb*Nb vector that stores the appropriately weighted sum of all previous RK steps for use in the final update step
        w: determines the weight that each K contributes to the final result
     
     */
    int Ns2 = int(pow(Ns,2));
    int Nb2 = int(pow(Nb,2));
    complex<double> im = {0,1};
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    for(int ks = 0; ks < Ns; ++ks){
                        for(int kb = 0; kb < Nb; ++kb){
                            for(int ls = 0; ls < Ns; ++ls){
                                for(int lb = 0; lb < Nb; ++lb){
                        
                                    int idx1 = IDX_8D(ks,ls,kb,lb,is,js,ib,jb);
                                    K[idx1] = -im*(h2G2_comm[idx1] + psi[idx1] + pi[idx1]);
                                    RK_result[idx1] += w*K[idx1];
                                }
							}
						}
					}
				}
			}
		}
	}
}


void G2_step_GW(vector<complex<double>> &h2G2_comm, vector<complex<double>> &psi,vector<complex<double>> &psi0, vector<complex<double>> &pi, vector<complex<double>> &K,vector<complex<double>> &RK_result, double w)
{
    int Ns2 = int(pow(Ns,2));
    int Nb2 = int(pow(Nb,2));
    complex<double> im = {0,1};
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    for(int ks = 0; ks < Ns; ++ks){
                        for(int kb = 0; kb < Nb; ++kb){
                            for(int ls = 0; ls < Ns; ++ls){
                                for(int lb = 0; lb < Nb; ++lb){
                        
                                    int idx1 = IDX_8D(ks,ls,kb,lb,is,js,ib,jb);
                                    K[idx1] = -im*(h2G2_comm[idx1] + psi[idx1] + pi[idx1] - psi0[idx1]);
                                    RK_result[idx1] += w*K[idx1];
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
