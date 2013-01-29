#include <iostream>
#include<cmath>
#include"lattice.h"
#include"atom.h"

#include <armadillo>


using namespace std;
using namespace arma;

int main()
{
    Lattice l(1,1,1,"Ar", 2);
    l.writeVMDfile("unitcell.xyz", "Argon");

    return 0;
}

