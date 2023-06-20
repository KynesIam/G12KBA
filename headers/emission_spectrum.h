//
//  emission_spectrum.h
//  
//
//  Created by Cian Reeves on 2/27/23.
//

#ifndef emission_spectrum_h
#define emission_spectrum_h
#include <vector>
#include <complex>
#include <fftw3.h>
#define PI 3.14159265359
extern int time_steps;
double gaussian(double t, double tp, double delta)
{
    return (1.0/(delta*sqrt(2.0*PI)))*exp(-pow((t-tp)/delta,2)/2);
}

void emission_spectrum(vector<complex<double>> &I, vector<vector<complex<double>>> &G_lss_t1t2, vector<vector<complex<double>>> &rho_tot, double delta, double tp, double w_max,double w_min)
{
    complex<double> im = {0,1};

    double w = w_min;
    complex<double> trG=0;
    complex<double> trG2=0;
    double dt = dt_fixed;
    double dw = .0125;
    while(w < w_max){
        cout<<w<<"\n";
        int idx=0;
        double t1 = 0;
        double t2 = 0;
        complex<double> temp=0;
        for(int i = 0; i < rho_tot.size(); ++i){
            t2=t1;
            for(int j = i; j < rho_tot.size(); ++j){
                    trG=0;
                    trG2=0;

                    for(int k = 0; k < Ns*Nb; ++k){
                        trG +=G_lss_t1t2[idx][k*(Ns*Nb+1)];
                    }
                    if(i==j){
                        
                        temp+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(trG);
                    }
                    else{
                        for(int k = 0; k < Ns*Nb; ++k){
                            trG2+=G_lss_t1t2[idx][k*(Ns*Nb+1)];
                        }
                        temp+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(trG*exp(-im*(t1-t2)*w) - conj(trG2)*exp(im*(t1-t2)*w));
                    }
                    idx+=1;
                t2+=dt;
            }
            t1+=dt;
        }
        I.push_back(temp);
        w += dw;
        cout<<double(w - w_min)/(w_max - w_min)*100<<"%\r";
        cout.flush();
    }
    cout<<w;
}

void full_spectrum_kspace(vector<vector<complex<double>>> &I, vector<vector<complex<double>>> &G_lss_t1t2,vector<vector<complex<double>>> &G_grt_t1t2, vector<vector<complex<double>>> &rho_tot, double delta, double tp, double w_max,double w_min)
{
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_2D(x,y)  (x*Ns+y)
    
    complex<double> im = {0,1};
    
    
    
    double dt = dt_fixed;
    double dw = .1;
    
    fftw_plan p, q;
    for(int ib = 0; ib < 2; ++ib)
    {
        cout<<"\n Computing emission spectrum  in (E,k) space for band "<<ib+1<<"\n";

        double w = w_min;
        int count = 0;

        while(w < w_max){
            int idx=0;
            double t1 = 0;
            double t2 = 0;
            vector<complex<double>> temp(Ns);

            for(int i = 0; i < rho_tot.size(); ++i){
                t2=t1;
                for(int j = i; j < rho_tot.size(); ++j){
                    vector<complex<double>> temp_lss(Ns),in_lss(Ns*Ns), out_lss(Ns*Ns),temp_grt(Ns),in_grt(Ns*Ns), out_grt(Ns*Ns);

                    p = fftw_plan_dft_2d(Ns,Ns, reinterpret_cast<fftw_complex*>(&in_lss[0]), reinterpret_cast<fftw_complex*>(&out_lss[0]), FFTW_FORWARD, FFTW_ESTIMATE);
                    for (int ks = 0; ks < Ns; ks++) {
                        for (int ls = 0; ls < Ns; ls++) {
                            
                            int idx1=IDX_4D(ls,ib,ks,ib);
                            int idx2=IDX_2D(ls,ks);
                            in_lss[idx2] = G_lss_t1t2[idx][idx1];
                        }
                    }
                    fftw_execute(p);
                    
                    p = fftw_plan_dft_2d(Ns,Ns, reinterpret_cast<fftw_complex*>(&in_grt[0]), reinterpret_cast<fftw_complex*>(&out_grt[0]), FFTW_FORWARD, FFTW_ESTIMATE);
                    for (int ks = 0; ks < Ns; ks++) {
                        for (int ls = 0; ls < Ns; ls++) {
                            
                            int idx1=IDX_4D(ls,ib,ks,ib);
                            int idx2=IDX_2D(ls,ks);
                            in_grt[idx2] = G_grt_t1t2[idx][idx1];
                        }
                    }
                    fftw_execute(p);
                    for (int k = 0; k < Ns; ++k) {
                        for (int l = 0; l < Ns; ++l) {
                            temp_grt[k] += out_grt[k*Ns+l];
                            temp_lss[k] += out_lss[k*Ns+l];

                        }
                    }
                    
                    if(i==j){

                        for(int k = 0; k < Ns; ++k){
                            
                            temp[k]+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(temp_lss[k]-temp_grt[k]);
                        }
                        
                    }
                    else{
                        for(int k = 0; k < Ns; ++k){
                            
                            temp[k]+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*((temp_lss[k]-temp_grt[k])*exp(-im*(t1-t2)*w) - conj(temp_lss[k] - temp_grt[k])*exp(im*(t1-t2)*w));
                        }
                    }
                    idx+=1;
                    t2+=dt;
                }
                t1+=dt;
            }
            
            if(ib == 1){
                for(int z = 0; z < temp.size(); ++z){
                    I[count][z]+=temp[z];

                }
            }
            else{
                I.push_back(temp);
            }
            ++count;

            w += dw;
            cout<<double(w - w_min)/(w_max - w_min)*100<<"%\r";
            cout.flush();
        }
    }
    fftw_destroy_plan(p);

}
void particle_spectrum_kspace(vector<vector<complex<double>>> &I, vector<vector<complex<double>>> &G_lss_t1t2,vector<vector<complex<double>>> &G_grt_t1t2, vector<vector<complex<double>>> &rho_tot, double delta, double tp, double w_max,double w_min)
{
#define IDX_4D(x,y,z,t)  ((x) + (y)*Ns + ((z) + (t)*Ns)*Nb*Ns)
#define IDX_2D(x,y)  (x*Ns+y)
    
    complex<double> im = {0,1};
    
    
    
    double dt = dt_fixed;
    double dw = .05;
    
    fftw_plan p, q;
    for(int ib = 0; ib < 2; ++ib)
    {
        cout<<"\n Computing emission spectrum  in (E,k) space for band "<<ib+1<<"\n";

        double w = w_min;
        int count = 0;

        while(w < w_max){
            int idx=0;
            double t1 = 0;
            double t2 = 0;
            vector<complex<double>> Iw(Ns);

            for(int i = 0; i < rho_tot.size(); ++i){
                t2=t1;
                for(int j = i; j < rho_tot.size(); ++j){
                    vector<complex<double>> Gkgrt(Ns),in_grt(Ns*Ns), out_grt(Ns*Ns);

                    
                    
                    p = fftw_plan_dft_2d(Ns,Ns, reinterpret_cast<fftw_complex*>(&in_grt[0]), reinterpret_cast<fftw_complex*>(&out_grt[0]), FFTW_FORWARD, FFTW_ESTIMATE);
                    for (int ks = 0; ks < Ns; ks++) {
                        for (int ls = 0; ls < Ns; ls++) {
                            
                            int idx1=IDX_4D(ls,ib,ks,ib);
                            int idx2=IDX_2D(ls,ks);
                            in_grt[idx2] = G_grt_t1t2[idx][idx1];
                        }
                    }
                    fftw_execute(p);
                    for (int k = 0; k < Ns; ++k) {
                        for (int l = 0; l < Ns; ++l) {
                            Gkgrt[k] += out_grt[k*Ns+l];
                        }
                    }
                    
                    if(i==j){
                        for(int k = 0; k < Ns; ++k){
                            
                            Iw[k]+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(-Gkgrt[k]);
                        }
                        
                    }
                    else{
                        for(int k = 0; k < Ns; ++k){
                            
                            Iw[k]+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(-(Gkgrt[k])*exp(-im*(t1-t2)*w) + conj(Gkgrt[k])*exp(im*(t1-t2)*w));
                        }
                    }
                    idx+=1;
                    t2+=dt;
                }
                t1+=dt;
            }
            
            if(ib == 1){
                for(int k = 0; k < Iw.size(); ++k){
                    I[count][k]+=Iw[k];

                }
            }
            else{
                I.push_back(Iw);
            }
            ++count;

            w += dw;
            cout<<double(w - w_min)/(w_max - w_min)*100<<"%\r";
            cout.flush();
        }
    }
    fftw_destroy_plan(p);

}


void emission_spectrum_full(vector<vector<complex<double>>> &I, vector<vector<complex<double>>> &G_lss_t1t2,double delta, double tp, double w_max,double w_min)
{
    complex<double> im = {0,1};

    double w = w_min;
    double dt = dt_fixed;
    double dw = .025;
    double eps=0;
    while(w < w_max){
        int idx=0;
        double t1 = 0;
        double t2 = 0;
        vector<complex<double>> temp(G_lss_t1t2[0].size());
        for(int i = 0; i < time_steps; ++i){
            t2=t1;
            for(int j = i; j < time_steps; ++j){
                  
                    if(i==j){
                        for(int k = 0; k < Ns*Nb*Ns*Nb; ++k){
                            temp[k]+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*G_lss_t1t2[idx][k];
                        }
                    }
                    else{
                        for(int k = 0; k < Ns*Nb*Ns*Nb; ++k){
                            
                            temp[k]+=dt*dt*gaussian(t1,tp,delta)*gaussian(t2,tp,delta)*(G_lss_t1t2[idx][k]*exp(-im*(t1-t2)*w ) - conj(G_lss_t1t2[idx][k])*exp(im*(t1-t2)*w));
                        }
                    }
                    idx+=1;
                t2+=dt;
            }
            t1+=dt;
        }
        I.push_back(temp);
        w += dw;
        cout<<double(w - w_min)/(w_max - w_min)*100<<"%\r";
        cout.flush();
    }
}
#endif /* emission_spectrum_h */
