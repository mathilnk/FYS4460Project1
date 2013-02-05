#include "verlet_solver.h"
using namespace std;

Verlet_solver::Verlet_solver(Lattice &molecule){
    this->molecule = molecule;
    current_time_step = 0;
}

void Verlet_solver::solve(double t_start,int timesteps, double dt, string filename){
    /*
      Filename should be without ending. If you want the files to be called test*.xyz, just write filename = test.
      This method solves for many timesteps, and for each timestep a xyz file is made.
      */

    double time = t_start;
    for(int i=0;i<timesteps;i++){
        time = t_start + dt*i;
        this->molecule = solve_one_time_step(time, dt, filename);
    }

}



Lattice Verlet_solver::solve_one_time_step(double t, double dt, string filename){
    /*
      Solves for one time step, and returns the Lattice at the next timestep.
      dt is the time_step_size,
      t is the time
      filename is the base for the filename the lattice is saved to (VMD style).
      if you want the filename to be test*.xyz, write string filename = test.
      */

    /////////////////////////////////////////////////////////////
    //gives the filename an ending, with the current_time_step//
    ///////////////////////////////////////////////////////////
    string current_time_step_string;
    stringstream out;
    out<<current_time_step;
    current_time_step_string = out.str();
    string filename_end = filename + current_time_step_string + ".xyz";
    //writes to file
    molecule.writeVMDfile(filename_end,"comment");

    ///////////////////////////////////////////////////////////////////////
    //calculates the new velocity and position with the Verlet algorithm//
    /////////////////////////////////////////////////////////////////////
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
        vec dr = (v + force(t)/(2*molecule.mass)*dt)*dt;
        cout<<dr<<endl;
        r_new = r + dr;
        a->position = r_new;
        a->velocity = v_new;
    }


    //update current_time_step;
    current_time_step++;
    return molecule;
}

vec Verlet_solver::force(double t){
    /*
      returns a vec of zeros
      */
    vec ze = zeros(3,1);
    return ze;
}



