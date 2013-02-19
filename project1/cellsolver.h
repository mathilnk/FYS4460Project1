#ifndef CELLSOLVER_H
#define CELLSOLVER_H
#include"lattice.h"
#include<math.h>
#include"cellcontainer.h"
#include"normal.hpp"
class CellSolver
{
public:
    CellSolver(int CellNx, int CellNy, int CellNz, int Nx, int Ny, int Nz, double b, double T, double r_cut, string element);
    //system parameters
    int CellNx;
    int CellNy;
    int CellNz;
    int CellN;
    int Nx;
    int Ny;
    int Nz;
    int N;
    double b;
    double T;
    double mass;
    double k;
    double sig, mean,pi, r_cut;
    double Lx,Ly,Lz;
    int seed,current_time_step;
    double volume;
    //Conversion factors
    double T0,L0;
    double dt;


    //initialize
    void initializeContainer();
    void makeLattice();
    void findPosAndMakeAtoms(vec posBase);
    //write to file:
    void writeToVMDfile(string filename, string comment, string element);
    void writeEnergyToFile(string filename);
    void writeMeasurementsToFile(string filename);
    //Solve_methods
    void solve_one_time_step(double t, double dt, string filename, bool writeVMD,bool writeMeasurements);
    void solve(double t_start, int timesteps, double dt, string filename, bool writeVMD=true,bool writeMeasurements=false);

    //forces
    void findForces();
    void cleanForces();
    vec force_between(vec r_1, vec r_2, bool sameCell=true);


    //update and measurements
    void updateCells();
    void updateMeanSquareDisplacement(vec r_new, vec r_old);
    double measureTemp();
    double measurePressure();

    //other
    double gauss(double s, double mean);
    vec findRadial();



    CellContainer* myContainer;
    string element;
    int counter;

    double d_energy;

    double pressure_sum;
    double d_pressure;
    double current_displacement;

    //results:
    vec kin_energy;
    vec pot_energy;
    vec temperature;
    vec pressure;
    vec displacement;
};

#endif // CELLSOLVER_H
