#include "thermostattest.h"

ThermostatTest::ThermostatTest()
{

    //runBandA_Always();
    runBUpDown(300,2,200,200);

}

void ThermostatTest::runBandA_Always(){
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
    string element = "Ar";
    double dt = 0.01;

    double numOfTimeSteps = 300;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string B_filename = "B_thermo";
    string A_filename = "A_thermo";
    bool Berendsen = true;
    bool Andersen = false;
    double T_bath = 100;
    CellSolver* bSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element);
    bSolver->thermostatON = true;
    CellSolver* aSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element);
    aSolver->thermostatON = true;

    bSolver->solve(0, numOfTimeSteps, dt, B_filename, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    Andersen = true;
    Berendsen = false;
    aSolver->solve(0, numOfTimeSteps, dt, A_filename, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
}

void ThermostatTest::runBUpDown(double first, double second, double first_timestep, double second_timestep){
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
    string element = "Ar";
    double dt = 0.01;

    double numOfTimeSteps;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string filename1 = "UpDown";
    //string filename2 = "UpDown2";
    //uble numOfTimeSteps;


    bool Berendsen = true;
    bool Andersen = false;
    double T_bath;
    int max_time_step = first_timestep+second_timestep;
    CellSolver* bSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,max_time_step, b, T, r_cut, element);
    bSolver->thermostatON = true;
    numOfTimeSteps = first_timestep;
    T_bath = first;
    bSolver->solve(0, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    cout<<"done wiht first"<<endl;
    T_bath = second;
    numOfTimeSteps = second_timestep;
    bSolver->solve(dt*first_timestep, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    cout<<"done with second"<<endl;
}
