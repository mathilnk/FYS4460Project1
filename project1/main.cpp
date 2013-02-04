#include <iostream>
#include<cmath>
#include"lattice.h"
#include"atom.h"

#include <armadillo>
#include"normal.hpp"
#include "verlet_solver.h"


using namespace std;
using namespace arma;

int main()
{
    int Ny,Nx,Nz;
    double T, mass, b;
    string element = "Ar";
    T = 100; //K
    mass =39.948*1.66e-27;
    b = 5.26; //Å
    //b=1;
    //T = 200;
    //mass = 1;
    Lattice l(2,2,2,"Ar", b, T,mass);
    //l.writeVMDfile("test.xyz", "comment");
    Verlet_solver v(l);
    char file[100] = "test0";
    //v.solve_one_time_step(0.1, 0.1, file);
    v.solve(0,10,0.1,file);



    return 0;
}

