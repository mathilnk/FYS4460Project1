#include "energytest.h"
using namespace arma;
EnergyTest::EnergyTest()
{
    double L0 = 3.405;
    double r_cut = 3*L0;
    int Nx = 8, Ny = 8, Nz = 8;
    int CellNx, CellNy, CellNz;
    double b = 5.26;
    double T = 100;
    CellNx = Nx*b/r_cut;
    //CellNx = 1;
    CellNy = CellNx;
    CellNz = CellNx;
    r_cut = Nx*b/(CellNx);
    double numOfTimeSteps = 1000;
//    CellNx = 1;
//    CellNy = 1;
//    CellNz = 1;
//    r_cut = Nx*b;



    //double maxTime = 1;
    string element = "Ar";
    string filenamebase = "test";
    string current_dt_string;
    stringstream out;
    bool writeVMD = true;
    bool writeMeasurements = true;
    double T_bath = 0;
    bool Andersen = false;
    bool Berendsen = false;



    //vec dt = linspace<vec> (0.1, 1, 10);
    vec dt = vec(3);
    dt<<0.006;
    cout<<dt<<endl;
    for(int i =0;i<dt.size(); i++){

        CellSolver* mySolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz, b, T, r_cut, element);


        cout<<"AT "<<i<<" experiment of "<<dt.size()<<" experiments"<<endl;
        //string filenamebase = "things";
        string current_dt_string;
        stringstream out;
        out<<i;
        current_dt_string = out.str();
        cout<<"her"<<endl;
        mySolver->solve(0,numOfTimeSteps, dt[i],filenamebase, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);  }


}
