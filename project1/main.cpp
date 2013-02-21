#include <iostream>
#include<cmath>
#include"lattice.h"
#include"atom.h"

#include <armadillo>
#include"normal.hpp"
#include "verlet_solver.h"
#include<time.h>
#include"cellsolver.h"
#include"atomnode.h"
#include"energytest.h"


using namespace std;
using namespace arma;

int main()
{
//    int Ny,Nx,Nz;
//    double T, mass, b;
//    string element = "Ar";
//    T = 100; //K
//    //mass =39.948*1.66e-27;
//    mass = 39.948; //amu
//    b = 5.26; //Å
//    //b=1;
//    //T = 200;
//    //mass = 1;
//    Lattice l(8,8,8,"Ar", b, T,mass);
//    //l.writeVMDfile("test.xyz", "comment");
//    Verlet_solver v(l);
//    string file = "ny_versjon";
//    //v.solve_one_time_step(0.1, 0.1, file);
//    clock_t start = clock();
//    v.solve(0,30,0.01,file);
//    clock_t stop = clock();
//    cout<<(stop-start)*1000/CLOCKS_PER_SEC<<" ms"<<endl;
//    int Nx,Ny,Nz;
//    Nx = 2;
//    Ny = 2;
//    Nz = 2;
//    double r_cut = 3;
////    CellContainer C = CellContainer(Nx,Ny,Nz,r_cut);
      vec pos =zeros(3,1);
      vec vel = zeros(3,1);
  //    pos<<5<<2<<5;
////    cout<<pos<<endl;
////    string e = "hei";
//      Atom * a = new Atom(pos,vel,"a");
//      Atom * b = new Atom(pos,vel,"b");
//    //C.putAtomInRightCell(a);
//    AtomNode* first = new AtomNode(a);
//    AtomNode * second = new AtomNode(new Atom(pos,vel, "2."));
//    first->addAtom(new Atom(pos,vel, "2."));
//    AtomNode* last = first->findLast();
////    last->nextAtom = new AtomNode(new Atom(pos,vel, "3."));
//      Cell test_cell = Cell();
//      test_cell.addAtom(a);
//      test_cell.addAtom(a);
//      test_cell.addAtom(a);
//      cout<<test_cell.numberOfAtoms<<endl;
//      test_cell.removeAtom(0);
//      cout<<test_cell.numberOfAtoms<<endl;
//      test_cell.removeAtom(0);
//      cout<<test_cell.numberOfAtoms<<endl;


      /*double L0 = 3.405;
      double r_cut = 3*L0;
      int Nx = 8, Ny = 8, Nz = 8;
      int CellNx, CellNy, CellNz;
      double b = 5.26;
      double T = 100;
      CellNx = Nx*b/r_cut;
      r_cut = Nx*b/(CellNx);


      string element = "Ar";*/

      //cells
      //CellSolver* mySolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz, b, T, r_cut, element);
      //string fileCell = "testCellForce";
      //string fileVerlet = "testVerlet";
      //mySolver->solve(0,300,0.006,fileCell);

      //double mass = 3;

      //Verlet
      int CellNx,Nx;
      double b, T,r_cut;
      string element = "Ar";
      Nx = 8;

      CellNx = 1;
      b = 5.28;
      T = 100;
      r_cut = b*Nx;
      double mass = 1;
      //Lattice mol = Lattice(Nx,Nx,Nx, element,b, T, mass);
      //Verlet_solver* vSolver = new Verlet_solver(mol);
      //string fileVerlet = "Verlet_solver_g";
      //vSolver->solve(0,1000,0.01,fileVerlet);
      //mol->writeVMDfile("testLattice.xyz", "kommentar");


              //      CellContainer* c = new CellContainer(8,8,8,2);
//      for(int i=0;i<12;i++){
//          pos<<i*0.5<<0<<i*0.3;
//          //Atom* A = new Atom(pos,vec,"Ar");
//          Atom * a = new Atom(pos,vel,"Ar");
//          c->addAtom(a);
//      }
      new EnergyTest();

//cout<<r_cut<<endl;

    //CellSolver *mySolver = new CellSolver(CellNx, CellNx, CellNx, Nx, Nx,Nx, b, T, r_cut,element);
    //mySolver->writeToVMDfile("exercise_a.xyz", "Exercise a. Visualize Argon solid",element);
    //mySolver->solve(0, 100,0.01,'Exercise_g', true, true);


//    AtomNode* first = new AtomNode(a);
//    string e2 = "YEY";
//    Atom * b = new Atom(pos,vel,e2);
//    first->nextAtom = new AtomNode(b);
//    cout<<first->isLast()<<endl;
//    AtomNode * c = first->findLast();
//    cout<<c->myAtom->element<<endl;
//    first->removeNodeAfterMe();
//    c = first->findLast();
//    cout<<c->myAtom->element<<endl;








    return 0;
}


