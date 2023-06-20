//
//  GKBA_collision_integral.h
//  
//
//  Created by Cian Reeves on 2/24/23.
//

#ifndef Self_energy_h
#define Self_energy_h
#include <complex>
#include <vector>
extern double U;

//void scnd_brn_SE(vector<vector<complex<double>>> &SE_lss_t1t2,vector<vector<complex<double>>> &SE_grt_t1t2,vector<vector<complex<double>>> &G_lss_t1t2,vector<vector<complex<double>>> &G_grt_t1t2,vector<complex<double>> &W){
//#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
//#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
//    for(int t_index = 0; t_index < G_lss_t1t2.size(); ++ t_index){
//        vector<complex<double>> SE_lss_temp(Ns*Ns*Nb*Nb);
//        vector<complex<double>> SE_grt_temp(Ns*Ns*Nb*Nb);
//
//        for(int is = 0; is < Ns; ++is){
//            for(int ib = 0; ib < Nb; ++ib){
//                for(int js = 0; js < Ns; ++js){
//                    for(int jb = 0; jb < Nb; ++jb){
//                        complex<double> C_grt = 0;
//                        complex<double> C_lss = 0;
//
//                        for(int rs = 0; rs < Ns; ++rs){
//                            for(int rb = 0; rb < Nb; ++rb){
//                                for(int ps = 0; ps < Ns; ++ps){
//                                    for(int pb = 0; pb < Nb; ++pb){
//                                        complex<double> B_grt = 0;
//                                        complex<double> B_lss = 0;
//
//                                        for(int ls = 0; ls < Ns; ++ls){
//                                            for(int lb = 0; lb < Nb; ++lb){
//                                                for(int qs = 0; qs < Ns; ++qs){
//                                                    for(int qb = 0; qb < Nb; ++qb){
//                                                        complex<double> A_grt = 0;
//                                                        complex<double> A_lss = 0;
//
//                                                        for(int ts = 0; ts < Ns; ++ts){
//                                                            for(int tb = 0; tb < Nb; ++tb){
//                                                                for(int ks = 0; ks < Ns; ++ks){
//                                                                    for(int kb = 0; kb < Nb; ++kb){
//                                                                        int idx2 = IDX_4D(ks,kb,ts,tb);
//                                                                        int idx3 = IDX_8D(ls,ps,lb,pb,is,ks,ib,kb);
//                                                                        int idx4 = IDX_8D(ts,js,tb,jb,qs,rs,qb,rb);
//
//                                                                        A_grt+=W[idx4]*W[idx3]*G_lss_t1t2[t_index][idx2];
//                                                                        A_lss+=W[idx4]*W[idx3]*G_grt_t1t2[t_index][idx2];
//                                                                    }
//                                                                }
//                                                            }
//                                                        }
//                                                        int idx5 = IDX_4D(ls,lb,qs,qb);
//                                                        B_grt+=-A_grt*conj(G_grt_t1t2[t_index][idx5]);
//                                                        B_lss+=-A_lss*conj(G_lss_t1t2[t_index][idx5]);
//                                                    }
//                                                }
//                                            }
//                                        }
//                                        int idx6 = IDX_4D(ps,pb,rs,rb);
//                                        C_grt+=-B_grt*conj(G_grt_t1t2[t_index][idx6]);
//                                        C_lss+=-B_lss*conj(G_lss_t1t2[t_index][idx6]);
//                                    }
//                                }
//                            }
//                        }
//                        int idx1 = IDX_4D(js,jb,is,ib);
//                        SE_grt_temp[idx1] = C_grt;
//                        SE_lss_temp[idx1] = C_lss;
//                    }
//                }
//            }
//        }
//        SE_grt_t1t2.push_back(SE_grt_temp);
//        SE_lss_t1t2.push_back(SE_lss_temp);
//    }
//}


void scnd_brn_SE(vector<vector<complex<double>>> &SE_lss_t2t1,vector<vector<complex<double>>> &SE_grt_t2t1,vector<vector<complex<double>>> &G_lss_t1t2,vector<vector<complex<double>>> &G_grt_t1t2,vector<complex<double>> &W){
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    for(int t_index = 0; t_index < G_lss_t1t2.size(); ++t_index){
        vector<complex<double>> SE_lss_temp(Ns*Ns*Nb*Nb);
        vector<complex<double>> SE_grt_temp(Ns*Ns*Nb*Nb);

        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        int idx1 = IDX_4D(js,jb,is,ib);
                        int idx2 = IDX_4D(is,ib,js,jb);
                        SE_grt_temp[idx1] = U*U*conj(G_grt_t1t2[t_index][idx2])*conj(G_grt_t1t2[t_index][idx2])*G_lss_t1t2[t_index][idx2];
                        SE_lss_temp[idx1] = U*U*conj(G_lss_t1t2[t_index][idx2])*conj(G_lss_t1t2[t_index][idx2])*G_grt_t1t2[t_index][idx2];
                    }
                }
            }
        }
        SE_grt_t2t1.push_back(SE_grt_temp);
        SE_lss_t2t1.push_back(SE_lss_temp);
    }
}


void scnd_brn_SE2(vector<vector<complex<double>>> &SE_lss_t2t1,vector<vector<complex<double>>> &SE_grt_t2t1,vector<vector<complex<double>>> &rho_lss_tot,vector<vector<complex<double>>> &Y,vector<complex<double>> &W){
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    for(int t_index = 0; t_index < rho_lss_tot.size(); ++t_index){
        vector<complex<double>> SE_lss_temp(Ns*Ns*Nb*Nb);
        vector<complex<double>> SE_grt_temp(Ns*Ns*Nb*Nb);
        vector<complex<double>> rho_grt(Ns*Ns*Nb*Nb);
        vector<complex<double>> Y_dag(Ns*Ns*Nb*Nb);

        vector<complex<double>> tmp1(Ns*Ns*Nb*Nb);
        vector<complex<double>> tmp2(Ns*Ns*Nb*Nb);
        vector<complex<double>> tmp3(Ns*Ns*Nb*Nb);
        vector<complex<double>> tmp4(Ns*Ns*Nb*Nb);
        
        for(int i = 0; i < Ns; ++i){
            for(int j = 0; j < Ns; ++j){
                Y_dag[i*Ns + j] = conj(Y[t_index][j*Ns + i]);
                if(i==j){
                    rho_grt[i*Ns+j] = 1.0-rho_lss_tot[t_index][i*Ns + j];
                }
                else{
                    rho_grt[i*Ns+j] = -rho_lss_tot[t_index][i*Ns+j];
                }
            }
        }
        matrix_mult(tmp1,rho_grt,Y[t_index]);
        matrix_mult(tmp2,Y_dag,rho_lss_tot[t_index]);
        
        matrix_mult(tmp3,rho_lss_tot[t_index],Y[t_index]);
        matrix_mult(tmp4,Y_dag,rho_grt);
        for(int is = 0; is < Ns; ++is){
            for(int ib = 0; ib < Nb; ++ib){
                for(int js = 0; js < Ns; ++js){
                    for(int jb = 0; jb < Nb; ++jb){
                        int idx1 = IDX_4D(js,jb,is,ib);
                        int idx2 = IDX_4D(is,ib,js,jb);
                        SE_lss_temp[idx1] = U*U*tmp1[idx2]*tmp2[idx1]*tmp2[idx1];
                        SE_grt_temp[idx1] = U*U*tmp3[idx2]*tmp4[idx1]*tmp4[idx1];

                    }
                }
            }
        }
        SE_grt_t2t1.push_back(SE_grt_temp);
        SE_lss_t2t1.push_back(SE_lss_temp);
    }
}



#endif
