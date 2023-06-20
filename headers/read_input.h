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

#ifndef read_input_h
#define read_input_h
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;


bool HF=false;
bool EHM=false;
bool sparse=true;
bool cluster=false;
double decay_rate=.7;
double t1=1;
double t2=0;
double t3=0;
double U=.5;
int time_steps=1000;
double dt_fixed=.02;
double quench_strength=1;
int Ns=2;
int Nq=0;
int Ns2=int(pow(Ns,2));
int Nb=1;
int Nb2=int(pow(Nb,2));
double a=1.0;
double lambda=2.0;
string SE="SB";
string q_type="none";
string interaction="full";
string decay_type="coulomb";
double coulomb_scaling=.5;
vector<double> epsilon{0,0};
vector<double> V;
vector<complex<double>> W;

vector<complex<double>> hsp;
vector<double> GW_evals;
double AS_rate=3;
double AS_midpoint=25;
double quench_rate=.2;
double quench_on_midpoint=50;
double quench_off_midpoint=55;
double t0=0;
double Tp=0;
double wp=0;
double E=0;
int dim=1;
double offset=0;
double alpha;
double G2_damping=0;
vector<vector<double>> r;
bool pb=false;
string read_var_val(string filename,string var)
{
    std::fstream myfile1;
    myfile1.open(filename);
    if(myfile1.fail()){
        cout<<"Missing input file 'input' \n";
    }
    int line_counter = 0;
    string line;
    while (getline(myfile1,line)){
        line_counter++;
        if (line.find(var) != string::npos){
            break;
        }
    }
    myfile1.close();
    
    std::fstream myfile2;
    myfile2.open(filename);
    string s;
    for(int i = 0; i < line_counter + 1; ++i){
        getline(myfile2,s);
    }
    if(s==""){
        cout<<"Missing value for" <<var<<"\n";
        return "0";
    }
    myfile2.close();
    return s;
}

int read_in_w(vector<complex<double>> &W,int Ns, int Nb)
{
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    string s;
    std::fstream myfile;
    myfile.open("W.txt");
    if(myfile.fail()){
        cout<<"Missing input file 'W.txt'\n";
        return 0;
    }
    for(int is = 0; is < Ns; ++is){
        for(int ib = 0; ib < Nb; ++ib){
            for(int js = 0; js < Ns; ++js){
                for(int jb = 0; jb < Nb; ++jb){
                    for(int ks = 0; ks < Ns; ++ks){
                        for(int kb = 0; kb < Nb; ++kb){
                            for(int ls = 0; ls < Ns; ++ls){
                                for(int lb = 0; lb < Nb; ++lb){
                                    int idx1 = IDX_8D(is,js,ib,jb,ls,ks,lb,kb);
                                    getline(myfile,s);

                                    if(interaction=="onsite"){
                                        if(js==is and ks==ls and js==ks){
                                            W[idx1] = stod(s);
                                        }
                                        else{
                                            W[idx1] = 0;
                                        }
                                    }
                                    else if(interaction=="extended"){
                                        if(js==is and ks==ls){
                                            W[idx1] = stod(s);
                                        }
                                        else{
                                            W[idx1] = 0;
                                        }
                                    }
                                    else if(interaction=="full"){
                                            W[idx1] = stod(s);
                                        
                                                                                
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    myfile.close();
    return 1;
}

int read_in_hsp(vector<complex<double>> &hsp, int Ns, int Nb)
{
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    string s;
    std::fstream myfile;
    myfile.open("hsp.txt");
    if(myfile.fail()){
        cout<<"Missing input file 'hsp.txt'\n";
        for(int i = 0; i < Ns*Ns*Nb*Nb; ++i){
            hsp.push_back(0);
        }
    }
    else{
        for(int i = 0; i < Ns*Ns*Nb*Nb; ++i){
            getline(myfile,s);
            hsp.push_back(stod(s));
        }
    }
    myfile.close();
    return 1;
}

int read_in_GW_evals(vector<double> &GW_evals,int N)
{
#define IDX_8D(x,y,z,t,m,n,q,w) ((x) + (y)*Ns + ((z) + (t)*Nb)*Ns2 + ((m)+(n)*Ns)*Ns2*Nb2 + ((q) + (w)*Nb)*Ns2*Ns2*Nb2)
    string s;
    std::fstream myfile;
    myfile.open("GW_evals.txt");
    if(myfile.fail()){
        cout<<"Missing input file 'GW_evals'\n";
        return 0;
    }
    for(int i = 0; i < N; ++i){
        getline(myfile, s);
        GW_evals.push_back(stod(s));
    }
    myfile.close();
    return 1;
    
}

void read_in_cfile(string filename,vector<complex<double>> &data)
{
    std::fstream myfile1;
    std::fstream myfile2;
    myfile1.open("imag_" + filename);
    myfile2.open("real_" + filename);
    double temp_i;
    double temp_r;
    for(int i = 0; i < data.size(); ++i){
        myfile1 >> temp_i;
        myfile2 >> temp_r;
        data[i] = {temp_r,temp_i};
    }
    myfile1.close();
    myfile2.close();
}

void read_in_cfile2D(string filename,vector<vector<complex<double>>> &data,int Ns2,int Nb2,int time_steps)
{
    std::fstream myfile1;
    std::fstream myfile2;
    myfile1.open("imag_" + filename);
    myfile2.open("real_" + filename);
    double temp_i;
    double temp_r;
    vector<complex<double>> temp_data(Ns2*Ns2*Nb2*Nb2);
    for(int j = 0; j < time_steps; ++j){
        for(int i = 0; i < temp_data.size(); ++i){
            
            myfile1 >> temp_i;
            myfile2 >> temp_r;
            temp_data[i] = {temp_r,temp_i};
            
        }
        data.push_back(temp_data);
    }
    myfile1.close();
    myfile2.close();
}

void assign_vals()
{
    
    t1=stod(read_var_val("input","t1"));
    t2=stod(read_var_val("input","t2"));
    t3=stod(read_var_val("input","t3"));
    U=stod(read_var_val("input","U"));
    decay_rate=stod(read_var_val("input","interaction_decay"));
    alpha=stod(read_var_val("input","hopping_decay"));
    coulomb_scaling=stod(read_var_val("input","coulomb_scaling"));

    E=stod(read_var_val("input","E_strength"));
    t0=stod(read_var_val("input","t0"));
    wp=stod(read_var_val("input","wp"));
    Tp=stod(read_var_val("input","Tp"));
    lambda=stod(read_var_val("input","wavelength"));
    quench_strength=stod(read_var_val("input","quench_strength"));
    Nq=stoi(read_var_val("input","Nq"));
    quench_on_midpoint=stod(read_var_val("input","quench_on_midpoint"));
    quench_off_midpoint=stod(read_var_val("input","quench_off_midpoint"));
    quench_rate=stod(read_var_val("input","quench_rate"));

    offset=stod(read_var_val("input","offset"));

    dim=stod(read_var_val("input","dim"));
    time_steps=stoi(read_var_val("input","time_steps"));
    dt_fixed=stod(read_var_val("input","dt"));
    Ns=stoi(read_var_val("input","Ns"));
    Nb=stoi(read_var_val("input","Nb"));
    Nb2=int(pow(Nb,2));
    Ns2=int(pow(Ns,2));


    AS_rate=stod(read_var_val("input","AS_rate"));
    AS_midpoint=stod(read_var_val("input","AS_midpoint"));
    
    G2_damping=stod(read_var_val("input","G2_damping"));

    
    if(read_var_val("input","HF") == "true" or read_var_val("input","HF") == "True"){
        HF=true;
    }
    else if(read_var_val("input","HF") == "false" or read_var_val("input","HF") == "False"){
        HF=false;
    }
    else{
        cout<<"Invalid specification of HF calculation parameter, performing full HF-GKBA calculation.\n";
        HF=false;
    }
    
    
    if(read_var_val("input","sparse") == "true" or read_var_val("input","sparse") == "True"){
        sparse=true;
    }
    else if(read_var_val("input","sparse") == "false" or read_var_val("input","sparse") == "False"){
        sparse=false;
    }
    else{
        cout<<"Invalid specification of parameter:sparse.  Defaulting to sparse matrix multiplication.\n";
        sparse=true;
    }
    
    
    if(read_var_val("input","EHM") == "true" or read_var_val("input","EHM") == "True"){
        EHM=true;
    }
    else if(read_var_val("input","EHM") == "false" or read_var_val("input","EHM") == "False"){
        EHM=false;
    }
    else{
        cout<<"Invalid specification for EHM parameter. performing full calculation for onsite model.\n";
        EHM=false;
    }
    
    if(read_var_val("input","cluster") == "true" or read_var_val("input","cluster") == "True"){
        cluster=true;
    }
    else if(read_var_val("input","cluster") == "false" or read_var_val("input","cluster") == "False"){
        cluster=false;
    }
    else{
        cout<<"Invalid specification for cluster parameter. performing calculation for NN hopping.\n";
        cluster=false;
    }
    
    if(read_var_val("input","pb") == "true" or read_var_val("input","pb") == "True"){
        pb=true;
    }
    else if(read_var_val("input","pb") == "false" or read_var_val("input","pb") == "False"){
        pb=false;
    }
    else{
        cout<<"Invalid specification for pb parameter. performing calculation with open boundary conditions.\n";
        pb=false;
    }
    
    if(read_var_val("input","interaction") == "extended" or read_var_val("input","interaction") == "Extended"){
        interaction = "extended";
    }
    else if(read_var_val("input","interaction") == "onsite" or read_var_val("input","interaction") == "Onsite"){
        interaction = "onsite";
    }
    else if(read_var_val("input","interaction") == "full" or read_var_val("input","interaction") == "Full"){
        interaction = "full";
    }
    else{
        cout<<"Invalid specification for interaction parameter. Using full interaction matrix.\n";
        interaction="full";
    }
        

    if(read_var_val("input","q_type")=="full"){
        q_type="full";
    }
    else if(read_var_val("input","q_type")=="pulse"){
        q_type="pulse";
    }
    else if(read_var_val("input","q_type")=="none"){
        q_type="none";
        quench_strength=0;
        Nq=0;
    }
    
    
    if(read_var_val("input","decay_type")=="exp"){
        decay_type="exp";
    }
    else if(read_var_val("input","decay_type")=="coulomb"){
        decay_type="coulomb";
    }
    else{
        cout<<read_var_val("input","decay_type")<<"\n";
        cout<<"Invalid specification of decay type, defaulting to exponential decay\n";
        decay_type="exp";
    }

    if(read_var_val("input","self_energy")=="SB" or read_var_val("input","self_energy")=="scnd_brn"){
        SE="scnd_brn";
    }
    else if(read_var_val("input","self_energy")=="GW"){
        SE="GW";
    }
    else{
        cout<<"Invalid specification of self energy, defaulting to second Born\n";
        SE="scnd_brn";
    }
    
    for(int i = 0; i < Nb2*Nb2*Ns2*Ns2; ++i){
        W.push_back(0);
    }
    for(int i = 0; i < Ns; ++i){
        r.push_back({(Ns/2) - i - .5,0,0});
    }
//    read_in_w(W,Ns,Nb);
//    read_in_hsp(hsp,Ns,Nb);
//  read_in_GW_evals(GW_evals,Ns*Nb);
}


#endif
