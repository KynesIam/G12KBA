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
#ifndef print_matrix_h
#define print_matrix_h
#include <vector>
#include <complex>
#include <iostream>
extern int Ns;
extern int Nb;
using namespace std;

void print_2idx_cmat(vector<complex<double>> &mat)
{
    for(int i = 0; i <Ns; ++i){
        for(int j = 0; j < Ns; ++j){
            for(int a = 0; a < Nb; ++ a){
                for(int b = 0;  b< Nb; ++b){
                    cout<<"["<<i<<","<<j<<","<<a<<","<<b<<"]:"<<mat[(i + a*Ns)*Nb*Ns + (j+b*Ns)]<<"\n";
                }
            }
        }
    }
    cout<<"\n";
}
void print_2idx_rmat(vector<double> &mat)
{
    for(int i = 0; i <Ns; ++i){
        for(int j = 0; j < Ns; ++j){
            for(int a = 0; a < Nb; ++ a){
                for(int b = 0;  b< Nb; ++b){
                    cout<<"["<<i<<","<<j<<","<<a<<","<<b<<"]:"<<mat[(i + a*Ns)*Nb*Ns + (j+b*Ns)]<<"\n";
                }
            }
        }
    }
    cout<<"\n";
}


void print_4idx_cmat(vector<complex<double>> &mat)
{
    int Ns2 = int(pow(Ns,2));
    int Nb2 = int(pow(Nb,2));
    for(int i = 0; i < Ns; ++i){
        for(int j = 0; j< Ns; ++j){
            for(int k = 0; k< Ns; ++k){
                for(int l = 0; l<Ns;++l){
                    for(int a = 0; a< Nb; ++a){
                        for(int b = 0; b< Nb; ++b){
                            for(int g = 0; g<Nb; ++g){
                                for(int d = 0; d<Nb;++d){
                                    cout<<"["<<i<<","<<j<<","<<k<<","<<l<<":"<<a<<","<<b<<","<<g<<","<<d<<"]:"<< mat[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2]<<"\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cout<<"\n";
}

void print_4idx_rmat(vector<double> &mat)
{
    int Ns2 = int(pow(Ns,2));
    int Nb2 = int(pow(Nb,2));
    for(int i = 0; i < Ns; ++i){
        for(int j = 0; j< Ns; ++j){
            for(int k = 0; k< Ns; ++k){
                for(int l = 0; l<Ns;++l){
                    for(int a = 0; a< Nb; ++a){
                        for(int b = 0; b< Nb; ++b){
                            for(int g = 0; g<Nb; ++g){
                                for(int d = 0; d<Nb;++d){
                                    cout<<"["<<i<<","<<j<<","<<k<<","<<l<<":"<<a<<","<<b<<","<<g<<","<<d<<"]:"<< mat[(i+j*Ns)*Ns2*Nb2 + k+Ns*l + (a + b*Nb)*Ns2*Ns2*Nb2 + (g + d*Nb)*Ns2]<<"\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cout<<"\n";
}

void print_rmat(vector<double> &mat)
{
    int N = int(sqrt(mat.size()));
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            cout<<mat[i*N+j]<<" ";
        }
        cout<<"\n";
    }
}

void print_cmat(vector<complex<double>> &mat)
{
    int N = int(sqrt(mat.size()));
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            cout<<mat[i*N+j]<<" ";
        }
        cout<<"\n";
    }
}


#endif
