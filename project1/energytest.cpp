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
    double T_bath = 1;
    bool Andersen = false;
    bool Berendsen = true;



    //vec dt = linspace<vec> (0.1, 1, 10);
    vec dt = vec(3);
    dt<<0.006;
    cout<<dt<<endl;
    vec T_vec(5);
    T_vec<<100;//<<50<<100<<300<<500;
    bool s = true;
    for(int i =0;i<1; i++){
        CellSolver* mySolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz, b, T, r_cut, element);
        T_bath = T_vec(i);
        T = T_bath;
        if(i==0){
            mySolver->slow = true;
        }
        cout<<"AT "<<i<<" experiment of "<<5<<" experiments"<<endl;
        //string filenamebase = "things";
        string current_string;
        stringstream out;
        out<<i;
        current_string = out.str();
        cout<<"her"<<endl;
        mySolver->solve(0,numOfTimeSteps, dt[i],filenamebase+current_string, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
        mySolver->slow = false;
    }

}
