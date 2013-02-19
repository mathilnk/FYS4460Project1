#include "energytest.h"
using namespace arma;
EnergyTest::EnergyTest()
{
    int CellNx = 12, CellNy = 12, CellNz = 12, Nx = 8, Ny=8, Nz=8;
    double b = 5.26;
    double T = 100;
    double r_cut = Nx*b/CellNx;
    double maxTime = 10;
    string element = "Ar";
    string filenamebase = "things";
    string current_dt_string;
    stringstream out;


    //vec dt = linspace<vec> (0.1, 1, 10);
    vec dt = vec(3);
    dt<<0.01;
    cout<<dt<<endl;
    for(int i =0;i<dt.size(); i++){
        CellSolver* mySolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz, b, T, r_cut, element);


        cout<<"AT "<<i<<" experiment of "<<dt.size()<<" experiments"<<endl;
        string filenamebase = "things";
        string current_dt_string;
        stringstream out;
        out<<i;
        current_dt_string = out.str();
        cout<<"her"<<endl;
        mySolver->solve(0,maxTime/dt[i], dt[i],filenamebase + current_dt_string, true, true);

    }


}
