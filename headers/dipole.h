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

#ifndef dipole_h
#define dipole_h
#include <iostream>
#include <vector>
#include <complex>
using namespace std;


void dipole(vector<double> &p,vector<vector<complex<double>>> &rho,vector<vector<double>> &r)
{
    /*
    Function to compute the dipole for a given set of coordinate vectors, for 3 dimensions the
    dipole is calculated along the z direction by default and for 1 dim along the one dimension available
    Args:
        p: vector to store the computed dipole, with length given by the total number of time steps
        rho: n_steps*(Ns*Ns*Nb*Nb) vector that holds the density matrix along the entire time evolution
        r: vector that stores coordinate vectors defining the lattice sites
    */
    int N = int(sqrt(rho[0].size()));
    for(int i = 0; i < rho.size(); ++i){
        for(int k = 0; k < N; ++k){
            if(dim == 3){
                p[i] += r[k][2]*rho[i][k*(N+1)].real();
            }
            if(dim == 1){
                p[i] += r[k][0]*rho[i][k*(N+1)].real();

            }
        }
    }
}
#endif
