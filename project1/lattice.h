#ifndef LATTICE_H
#define LATTICE_H
#include<armadillo>
#include<facecube.h>
#include<vector>
#include<iostream>
#include<fstream>
#include"atom.h"
using namespace std;
using namespace arma;
class Lattice
{
public:
    Lattice(int Nx, int Ny, int Nz, string element, double b);
    vector<Atom*> allAtoms;
    void writeVMDfile(const char* Filename, string comment);
    void findPosAndMakeAtoms(vec posBase);

private:
    int Nx,Ny,Nz,N;
    string element;
    double b;
    void makeLattice();
    int numberOfAtoms;



};

#endif // LATTICE_H
