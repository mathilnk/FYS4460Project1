#include "verlet_solver.h"
using namespace std;

Verlet_solver::Verlet_solver(Lattice &molecule){
    this->molecule = molecule;
    current_time_step = 0;
}

void Verlet_solver::solve(double t_start,int timesteps, double dt, char* filename){
    /*
      Filename should be a large char *, with room for num.xyz
      This method solves for many timesteps, and for each timestep a xyz file is made.
      */
    char new_filename[100];
    double time = t_start;
    for(int i=0;i<timesteps;i++){
        time = t_start + dt*i;
        strcpy(new_filename, filename);
        this->molecule = solve_one_time_step(time, dt, new_filename);
    }

}



Lattice Verlet_solver::solve_one_time_step(double t, double dt, char* filename){
    char buffer[100];
    char time_step_char[64];
    char end_xyz[46] = ".xyz";
    sprintf(time_step_char, "%d", current_time_step);
    strcat(filename, time_step_char);
    strcat(filename, end_xyz);
    molecule.writeVMDfile(filename,"comment");
    int length = molecule.numberOfAtoms;
    vec v_new(3);
    vec v(3);
    vec r_new(3);
    vec r(3);
    int n=3;
    for(int i=0;i<length;i++){
        Atom* a = molecule.allAtoms[i];
        v = a->velocity;
        v_new = v + force(t)/(2*molecule.mass)*dt + force(t+dt)/(2*molecule.mass)*dt;
        r = a->position;
        r_new = r + (v + force(t)/(2*molecule.mass)*dt)*dt;
        a->position = r;
        a->velocity = v;
    }


    //writeVMDfile_Verlet(filename, 'comment');
    current_time_step++;
    return molecule;
}

vec Verlet_solver::force(double t){
    vec ze = zeros(3,1);
    return ze;
}


void Verlet_solver::writeVMDfile_Verlet(char* filename, string comment){
//    char buffer[100];
//    char time_step_char[64];
//    char end_xyz[46] = ".xyz";
//    sprintf(time_step_char, "%d", current_time_step);
//    strcat(filename, time_step_char);
//    strcat(filename, end_xyz);
//    molecule.writeVMDfile(filename,comment);
//    current_time_step++;
}
