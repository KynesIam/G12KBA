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

#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/read_input.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/emission_spectrum.h"

#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/RK_solver.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/init_G.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/AS.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/dipole.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/write_to_file.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/print_matrix.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/two_body.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/Off_diagonals.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/Self_Energy.h"
#include "/Users/cianreeves/Desktop/Code/local copy G12KBA1.0/src/headers/HF.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>

using namespace std;

extern int Ns;
extern int Nb;
extern int Ns2;
extern int Nb2;
//extern vector<complex<double>> W;
//extern vector<double> GW_evals;
//extern vector<complex<double>> hsp;
extern bool HF;
extern int time_steps;
extern vector<vector<double>> r;
int main()
{
    assign_vals();
    
    vector<complex<double>> G2xyxy(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> G2xxxx(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> rho(Ns*Ns*Nb*Nb);
    vector<vector<complex<double>>> rho_diag;
    vector<vector<complex<double>>> rho_diag_AS;
    
    vector<complex<double>> W(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> h(Ns*Ns*Nb*Nb);


    h_sp_Hubbard_nd(h,r,  t1);
    cout<<"Single particle matrix at equilibrium:\n";
    print_cmat(h);
    make_Hubbard_w(W,U,r);

    rho0(rho,Ns*Nb/2,h);
    cout<<"Non-interacting equilibrium Green's function:\n";

    print_cmat(rho);
    if(U != 0){
        Adiabatic_switching(rho,G2xxxx,G2xyxy,W,h,rho_diag_AS);
        cout<<"Equilibrium Green's function:\n";
        print_cmat(rho);
    }
    
    
    RK4(rho,G2xxxx,G2xyxy,rho_diag,W,h,r);

    vector<complex<double>> I;

    vector<vector<complex<double>>> G_lss_t1t2;
    vector<vector<complex<double>>> G_grt_t1t2;
//
    G_less_t1t2_full(G_lss_t1t2,G_grt_t1t2, rho_diag,W,h);
//    particle_spectrum_kspace(I,G_lss_t1t2,G_grt_t1t2,rho_diag,10.0,time_steps*dt_fixed/2.0,10,-10);
//    particle_spectrum_kspace(I,G_lss_t1t2,G_grt_t1t2,rho_diag,10.0,,6.5,-6.5);

    emission_spectrum(I,G_lss_t1t2,rho_diag,12.0,time_steps*dt_fixed/2.0,2.5,-2.5);
    write_to_cfile("full_emission.txt",I);
    
    vector<double> p(rho_diag.size());
    dipole(p,rho_diag,r);

    write_to_cfile2D("density_matrix.txt",rho_diag);
    write_to_rfile("dipole.txt",p);

    return 0;
}

