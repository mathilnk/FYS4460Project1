#ifndef VERLET_SOLVER_H
#define VERLET_SOLVER_H
#include"lattice.h"

class Verlet_solver
{
public:
    Verlet_solver(Lattice &molecule);
    Lattice molecule;
    Lattice solve_one_time_step(double t, double dt, string filename);
    vec force(double t);
    void solve(double t_start,int timesteps, double dt, string filename);
private:
    int current_time_step;
    void writeVMDfile_Verlet(char*filename, string comment);
};

#endif // VERLET_SOLVER_H
