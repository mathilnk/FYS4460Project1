#ifndef LATTICE_H
#define LATTICE_H
#include<armadillo>
#include<facecube.h>
#include<vector>
#include<iostream>
#include<fstream>
#include"atom.h"
#include<random>
using namespace std;
using namespace arma;
class Lattice
{
public:
    Lattice(int Nx, int Ny, int Nz, string element, double b, double T);
    vector<Atom*> allAtoms;
    void writeVMDfile(const char* Filename, string comment);
    void findPosAndMakeAtoms(vec posBase);
    void findPosAndMakeAtoms2(vec posBase);
    void makeEndAtoms();
    double T;

private:
    int Nx,Ny,Nz,N;
    string element;
    double b;
    void makeLattice();
    int numberOfAtoms;
    double gauss(double s, double m);
    double s,m,k;
    double pi;



};

#endif // LATTICE_H
