#include <iostream>
#include<cmath>
#include"lattice.h"
#include"atom.h"

#include <armadillo>


using namespace std;
using namespace arma;

int main()
{

    //vec a = zeros(4);
    //vec p = zeros(2);
    //p[0] = 1;
    //p[1] = 2;
    //cout << a << endl;
    //a[2] = 3;
    //cout << a << endl;
    //mat b = zeros(5,5);

    //cout << b << endl;
    //b.save("test.mat", raw_ascii);
    //mat c;
    ///c.load("test2.mat");

    //cout << c << endl;
  string str = "Ar";
  int b=2;
    cout<<b*mat("0 0 0 0 0 0;0.5 0.5 0 0 0 0;0 0.5 0.5 0 0 0;0.5 0 0.5 0 0 0;")<<endl;

Lattice l(1,1,1,str, 2);
l.get_vmd("Halla");



    return 0;
}

