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


#ifndef write_to_file_h
#define write_to_file_h

#include <fstream>
#include <vector>
#include <complex>
using namespace std;

void write_to_rfile(string filename, vector<double> data)
{
    std::ofstream myfile;
    myfile.open(filename);
    for(int i = 0; i<(data.size()); i++){
            myfile << data[i] << '\n';
    }
    myfile.close();
}

void write_to_cfile(string filename,vector<complex<double> > data)
{
    std::ofstream myfile1;
    std::ofstream myfile2;
    myfile1.open("imag_" + filename);
    myfile2.open("real_" + filename);
    for(int i = 0; i<(data.size()); i++){
        myfile1 <<data[i].imag()<<" ";
        myfile2 <<data[i].real()<<" ";
    }
    myfile1<<"\n";
    myfile2<<"\n";

    myfile1.close();
    myfile2.close();
}

void write_to_cfile2D(string filename,vector<vector<complex<double> > > data)
{
    std::ofstream myfile1;
    std::ofstream myfile2;
    myfile1.open("imag_" + filename);
    myfile2.open("real_" + filename);

	for(int j = 0; j < data[0].size(); ++j){
        for(int i = 0; i<(data.size()); i++){
	        myfile1 <<data[i][j].imag()<<" ";
	        myfile2 <<data[i][j].real()<<" ";
        }
        myfile1<<"\n";
        myfile2<<"\n";
	}
    myfile1.close();
    myfile2.close();
}




#endif
