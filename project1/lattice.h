#ifndef LATTICE_H
#define LATTICE_H
#include<armadillo>
#include<facecube.h>
using namespace std;
using namespace arma;
class Lattice
{
public:
    Lattice(int Nx, int Ny, int Nz, string element, double b);
    void get_vmd(string comment);
    vector<FaceCube>* std_lattice;
    mat arma_lattice;
    string element;
    double b;

private:
    int Nx,Ny,Nz,N;
    void makeLattice();

};

#endif // LATTICE_H
