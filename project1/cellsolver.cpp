#include "cellsolver.h"



using namespace std;

CellSolver::CellSolver(int CellNx, int CellNy, int CellNz, int Nx, int Ny, int Nz, double b, double T, double r_cut, string element){
    this->CellNx = CellNx;
    this->CellNy = CellNy;
    this->CellNz = CellNz;
    this->CellN = CellNx*CellNy*CellNz;
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;

    T0 = 119.74; //Kelvin
    L0 = 3.405; //AAngstrom
    this->b = b/L0;
    this->T = T;
    this->mass = 1;
    k = 8.617e-5;//ev/K
    this->pi = 3.1415;
    this->seed = 123456754;
    this->Lx = (Nx)*this->b;
    this->Ly = (Ny)*this->b;
    this->Lz = (Nz)*this->b;
    this->r_cut = r_cut/L0;
    this->element = element;
    this->mean = 0;
    this->sig= sqrt(T/T0);
    initializeContainer();
    current_time_step = 0;
    counter = 0;
    volume = Nx*Ny*Nz*this->b*this->b*this->b;
    d_pressure = 0;
    current_displacement=0;
    //displacement = 0;





}




void CellSolver::solve(double t_start, int timesteps, double dt, string filename, bool writeVMD, bool writeMeasurements){
    /*
      Filename should be without ending. If you want the files to be called test*.xyz, just write filename = test.
      This method solves for many timesteps, and for each timestep a xyz file is made.
      */

    double time = t_start;

    this->kin_energy = zeros(timesteps,1);
    this->pot_energy = zeros(timesteps,1);
    this->pressure = zeros(timesteps,1);
    this->temperature = zeros(timesteps,1);
    this->displacement = zeros(timesteps,1);
    //displacement = 0;


    ///////////////////////////////////////////////////
    //before the first timestep: calculate the forces//
    ///////////////////////////////////////////////////
    clock_t start = clock();

    findForces();


    for(int i=0;i<timesteps;i++){
        cout<<"tidssteg: "<<current_time_step<<endl;
        time = t_start + dt*i;
        kin_energy[i]=0;
        pot_energy[i]=0;
        d_energy=0;
        d_pressure =0;

        solve_one_time_step(time, dt, filename, writeVMD, writeMeasurements);

    }
    clock_t stop = clock();

    if(writeMeasurements){
        cout<<filename<<endl;
        writeMeasurementsToFile(filename+"_energy.txt");

    }
    cout<<(stop-start)/CLOCKS_PER_SEC<< " s for Cell solver"<<endl;
}

void CellSolver::solve_one_time_step(double t, double dt, string filename, bool writeVMD, bool writeMeasurements){
    /*
      Solves for one time step, and returns nothing.
      dt is the time_step_size,
      t is the time
      filename is the base for the filename the lattice is saved to (VMD style).
      if you want the filename to be test*.xyz, write string filename = test.
      */

    /////////////////////////////////////////////////////////////
    //gives the filename an ending, with the current_time_step//
    ///////////////////////////////////////////////////////////
    this->dt = dt;
    string current_time_step_string;
    stringstream out;
    out<<current_time_step;
    current_time_step_string = out.str();
    string filename_end = filename + current_time_step_string + ".xyz";
    //finds the radial distribution function
    vec g = findRadial();
    cout<<g<<endl;
    //writes to file
    if(writeVMD){
        myContainer->writeVMDfile(filename_end,"comment", element);
    }
    //cout<<"0 "<<d_energy<<endl;

    ///////////////////////////////////////////////////////////////////////
    //calculates the new velocity and position with the Verlet algorithm//
    /////////////////////////////////////////////////////////////////////

    int n_atoms = myContainer->numberOfAtoms;
    vec v_new(3);
    vec v(3);
    vec r_new(3);
    vec r(3);
    vec dr(3);
    vec v_half(3);
    double dx,dy,dz;
    double Lx = Nx*b;
    double Ly = Ny*b;
    double Lz = Nz*b;
    Cell currentCell;
    AtomNode* currentNode;
    Atom* atm;
    pressure_sum = 0;

    //////////////////////
    //loop to find r_new//
    //////////////////////


    for(int j=0;j<CellN;j++){

        //Atom* atm = molecule.allAtoms[i];
        currentCell = myContainer->myCells[j];
        for(int i=0;i<currentCell.numberOfAtoms;i++){
            atm = currentCell.myAtoms[i];
            v = atm->velocity;
            r = atm->position;
            //v_half = v + force_on(r, i)/(2)*dt;
            v_half = v + atm->force/2*dt;
            r_new = r + v_half*dt;
            v_new = v_half; //saves v_half so I dont have to calculate it again in the loop to find v_new

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Periodic boundary conditions (assume that a particle wont go further than 1000 systems away in the negative direction//
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            updateMeanSquareDisplacement(r_new,r);
            r_new[0] = fmod((r_new[0]+1000*Lx), Lx);
            r_new[1] = fmod((r_new[1]+1000*Ly), Ly);
            r_new[2] = fmod((r_new[2]+1000*Lz), Lz);



            atm->position = r_new;
            atm->velocity = v_new;
        }
    }



    ///////////////////////////////////////////
    //finds the forces with the new positions//
    ///////////////////////////////////////////
    //cout<<"before forces"<<endl;
    updateCells();
    findForces();
    //cout<<"before 2 loop"<<endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //loop to find v_new, now the atoms all have a new position. v_half is saved in atm->velocity//
    ///////////////////////////////////////////////////////////////////////////////////////////////
    for(int j=0;j<CellN;j++){
        currentCell = myContainer->myCells[j];
        for(int i=0;i<currentCell.numberOfAtoms;i++){
            atm = currentCell.myAtoms[i];
            v_half = atm->velocity;
            v_new = v_half + atm->force/2*dt;
            atm->velocity = v_new;
            d_energy += 0.5*(pow(atm->velocity[0],2) + pow(atm->velocity[1],2) + pow(atm->velocity[2],2));
            kin_energy[current_time_step]+=0.5*(pow(atm->velocity[0],2) + pow(atm->velocity[1],2) + pow(atm->velocity[2],2));
        }
    }




    double temp = measureTemp();
    d_pressure = myContainer->numberOfAtoms*temp/volume + 1/(6*volume)*pressure_sum;
    cout<<"Temperature: "<<temp*T0<<endl;
    cout<<"Pressure: "<<measurePressure()<<endl;
    cout<<"Kinetic energy: "<<kin_energy[current_time_step]<<endl;
    cout<<"Potential energy: "<<pot_energy[current_time_step]<<endl;

    cout<<"displacement: "<<displacement[current_time_step]<<endl;
    cout<<myContainer->numberOfAtoms<<endl;
    current_time_step++;//update current_time_step;

    return;
}


void CellSolver::updateCells(){
    Cell currentCell;
    CellContainer* buffer = new CellContainer(CellNx, CellNy, CellNz, r_cut);

    for(int j=0;j<CellN;j++){
        currentCell = myContainer->myCells[j];
        for(int i=0;i<currentCell.numberOfAtoms;i++){
            buffer->addAtom(currentCell.myAtoms[i]);
        }
        currentCell.myAtoms.clear();
    }
    myContainer = buffer;

}

void CellSolver::cleanForces(){

    Cell currentCell;
    AtomNode* currentNode;
    Atom* atm;
    vec z = zeros(3,1);
   for(int j=0;j<CellN;j++){
        currentCell = myContainer->myCells[j];
        for(int i=0;i<currentCell.numberOfAtoms;i++){
            atm = currentCell.myAtoms[i];
            atm->force = z;
        }
    }
}





void CellSolver::findForces(){
    cleanForces();
    Cell currentCell;
    AtomNode* currentNode1, *currentNode2;
    Atom* atom_1, *atom_2;
    ////////////////////////////////////////////////////////////////
    //first find the forces between the particles in the same cell//
    ////////////////////////////////////////////////////////////////
    for(int j=0;j<CellN;j++){
        currentCell = myContainer->myCells[j];

        for(int i=0;i<currentCell.numberOfAtoms;i++){
            atom_1 = currentCell.myAtoms[i];

            for(int k=i+1;k<currentCell.numberOfAtoms;k++){
                atom_2 = currentCell.myAtoms[k];
                vec f = force_between(atom_1->position,atom_2->position);
                atom_1->force = atom_1->force + f;
                atom_2->force = atom_2->force - f;
            }
        }
    }








    ///////////////////////////////////////////////////////////////////////
    //Now I have to calculate the forces between atoms in different cells//
    ///////////////////////////////////////////////////////////////////////
    Cell currentCell1, currentCell2;
    Col<int> myNeighbors;

    for(int j=0;j<CellN;j++){
        myNeighbors = myContainer->findMyNeighbors(j);
        currentCell1 = myContainer->myCells[j];
        for(int i=0;i<26;i++){
            currentCell2 = myContainer->myCells[myNeighbors[i]];

            for(int k=0;k<currentCell1.numberOfAtoms;k++){
                atom_1 = currentCell1.myAtoms[k];
                for(int p=0;p<currentCell2.numberOfAtoms;p++){
                    atom_2 = currentCell2.myAtoms[p];
                    vec f = force_between(atom_1->position,atom_2->position, false);
                    atom_1->force = atom_1->force + f;
                    //atom_2->force = atom_2->force - f;

                }
            }
        }
    }
}



vec CellSolver::force_between(vec r_1, vec r_2, bool sameCell){
    /*
      finds the force on r_1 from r_2 and vica versa. the force on r_2 is the same vector, but with different sign.Dimensionless
      */
    vec r = r_1-r_2;
    double Lx = Nx*b;
    double Ly = Ny*b;
    double Lz = Nz*b;
    double new_r_x_min = r[0] - Lx;
    double new_r_y_min = -(Ly - r[1]);
    double new_r_z_min = -(Lz - r[2]);
    double new_r_x_plus = r[0] + Lx;
    double new_r_y_plus = (Ly + r[1]);
    double new_r_z_plus = (Lz + r[2]);




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Complicated expressions. They choose the components for r(distance between particles) by choosing the smallest of (x_i - x_j + dL) where dL can be (-L,0,L)//
    //Called minimum image convention. (we have periodic boundary conditions.                                                                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    r[0] = (fabs(new_r_x_min)<fabs(new_r_x_plus))*(fabs(new_r_x_min)<fabs(r[0]))*new_r_x_min  +  (fabs(new_r_x_plus)<fabs(new_r_x_min))*(fabs(new_r_x_plus)<fabs(r[0]))*new_r_x_plus    +   (fabs(r[0])<fabs(new_r_x_plus))*(fabs(r[0])<fabs(new_r_x_min))*r[0];
    r[1] = (fabs(new_r_y_min)<fabs(new_r_y_plus))*(fabs(new_r_y_min)<fabs(r[1]))*new_r_y_min  +  (fabs(new_r_y_plus)<fabs(new_r_y_min))*(fabs(new_r_y_plus)<fabs(r[1]))*new_r_y_plus    +   (fabs(r[1])<fabs(new_r_y_plus))*(fabs(r[1])<fabs(new_r_y_min))*r[1];
    r[2] = (fabs(new_r_z_min)<fabs(new_r_z_plus))*(fabs(new_r_z_min)<fabs(r[2]))*new_r_z_min  +  (fabs(new_r_z_plus)<fabs(new_r_z_min))*(fabs(new_r_z_plus)<fabs(r[2]))*new_r_z_plus    +   (fabs(r[2])<fabs(new_r_z_plus))*(fabs(r[2])<fabs(new_r_z_min))*r[2];


    double length = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    length = std::max(0.9, length);
    vec f = -24*(1/pow(length,4)-2/pow(length,7))*r;
    if(sameCell){
        pot_energy[current_time_step] -=4*(1/pow(length,6) - 1/pow(length,3));
    }else{
        pot_energy[current_time_step] -= 2*(1/pow(length,6) - 1/pow(length,3));
    }
    d_energy-=4*(1/pow(length,6) - 1/pow(length,3));
    //cout<<4*(1/pow(length,6) - 1/pow(length,3))<<endl;
    pressure_sum = pressure_sum + f[0]*r[0] + f[1]*r[1] + f[2]*r[2];
    return f;
}



void CellSolver::initializeContainer(){
    myContainer = new CellContainer(CellNx, CellNy, CellNz, r_cut);
    makeLattice();
}


void CellSolver::makeLattice(){

    vec posBase(3);
    for(int k=0;k<Nz;k++){
        for(int j=0;j<Ny;j++){
            for(int i=0;i<Nx; i++){

                posBase<<i*b<<j*b<<k*b;
                findPosAndMakeAtoms(posBase);

            }
        }
    }
}


void CellSolver::writeToVMDfile(string filename, string comment, string element){
    myContainer->writeVMDfile(filename,comment,element);
}

void CellSolver::findPosAndMakeAtoms(vec posBase){
    /*
    given a position base x,y,z, this method creates four atoms as described in the exercise.
    The result is a face-centered cubic lattice
    */
    vec posxy(3);
    vec posyz(3);
    vec poszx(3);
    vec velxy(3);
    vec velyz(3);
    vec velzx(3);
    vec velbase(3);

    posxy<<0.5*b<<0.5*b<<0;
    posyz<<0<<0.5*b<<0.5*b;
    poszx<<0.5*b<<0<<0.5*b;
    posxy+=posBase;
    posyz+=posBase;
    poszx+=posBase;

    double s = sig;
    velxy<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    velyz<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    velzx<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    velbase<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    vec velo(3);
    velo<<0<<0<<0;
    Atom *xy = new Atom(posxy,velxy, element);
    Atom *yz = new Atom(posyz,velyz, element);
    Atom *zx = new Atom(poszx,velzx, element);
    Atom *base = new Atom(posBase,velbase, element);


    myContainer->addAtom(xy);

    myContainer->addAtom(yz);

    myContainer->addAtom(zx);

    myContainer->addAtom(base);



}

double CellSolver::measureTemp(){
    double temp = 2*kin_energy[current_time_step]/(3*myContainer->numberOfAtoms);
    temperature[current_time_step] = temp;
    return temp;
}

double CellSolver::measurePressure(){
    /*
      should rune measureTemp first
      */
    pressure[current_time_step] =  myContainer->numberOfAtoms*temperature[current_time_step]/volume + 1/(6*volume)*pressure_sum;
    return pressure[current_time_step];
}



void CellSolver::updateMeanSquareDisplacement(vec r_new, vec r_old){
    current_displacement += (r_new[0] - r_old[0])*(r_new[0] - r_old[0]) + (r_new[1] - r_old[1])*(r_new[1] - r_old[1]) + (r_new[2] - r_old[2])*(r_new[2] - r_old[2])/myContainer->numberOfAtoms;
    displacement[current_time_step] = current_displacement;
}
void CellSolver::writeMeasurementsToFile(string filename){
    ofstream file;
    string newfilename = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460Project1/results/" + filename;
    const char* filename_ch = newfilename.c_str();
    cout<<filename_ch<<endl;
    file.open(filename_ch);
    file<<"dt: "<<dt<<" CellNx: "<<CellNx<<" b: "<<b*L0<<" noParticles: "<<myContainer->numberOfAtoms<<endl;
    file<<"Kinetic energy     Potential energy   Total energy   Temperature         Pressure   Average displacement"<<endl;
    for(int i=0;i<current_time_step;i++){
        file<<kin_energy[i]<<"            "<<pot_energy[i]<<"               "<<kin_energy[i] + pot_energy[i]<<"         "<<temperature[i]<<"         "<<pressure[i]<<"         "<<displacement[i]<<endl;
    }
}

void CellSolver::writeEnergyToFile(string filename){
    ofstream file;
    const char* filename_ch = filename.c_str();
    cout<<filename_ch<<endl;
    file.open(filename_ch);
    for(int i=0;i<current_time_step;i++){
        file<<kin_energy[i]<<"      "<<pot_energy[i]<<"     "<<kin_energy[i] + pot_energy[i]<<endl;
    }

}

vec CellSolver::findRadial(){
    /*
      finds the radial distribution g(r)
      */
    double numOfBins = CellNx*10;
    vec bins = zeros(numOfBins, 1);
    Atom * atm1;
    Atom * atm2;
    double distance;
    double maxlength = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    int index;
    double binSize = maxlength/numOfBins;
    int maxIndex=0;
    vec r;
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        atm1 = myContainer->allAtoms[i];
        for(int j=i;j<myContainer->numberOfAtoms;j++){
            atm2 = myContainer->allAtoms[j];
            r = atm1->position - atm2->position;

            double new_r_x_min = r[0] - Lx;
            double new_r_y_min = -(Ly - r[1]);
            double new_r_z_min = -(Lz - r[2]);
            double new_r_x_plus = r[0] + Lx;
            double new_r_y_plus = (Ly + r[1]);
            double new_r_z_plus = (Lz + r[2]);




            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Complicated expressions. They choose the components for r(distance between particles) by choosing the smallest of (x_i - x_j + dL) where dL can be (-L,0,L)//
            //Called minimum image convention. (we have periodic boundary conditions.                                                                                    //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            r[0] = (fabs(new_r_x_min)<fabs(new_r_x_plus))*(fabs(new_r_x_min)<fabs(r[0]))*new_r_x_min  +  (fabs(new_r_x_plus)<fabs(new_r_x_min))*(fabs(new_r_x_plus)<fabs(r[0]))*new_r_x_plus    +   (fabs(r[0])<fabs(new_r_x_plus))*(fabs(r[0])<fabs(new_r_x_min))*r[0];
            r[1] = (fabs(new_r_y_min)<fabs(new_r_y_plus))*(fabs(new_r_y_min)<fabs(r[1]))*new_r_y_min  +  (fabs(new_r_y_plus)<fabs(new_r_y_min))*(fabs(new_r_y_plus)<fabs(r[1]))*new_r_y_plus    +   (fabs(r[1])<fabs(new_r_y_plus))*(fabs(r[1])<fabs(new_r_y_min))*r[1];
            r[2] = (fabs(new_r_z_min)<fabs(new_r_z_plus))*(fabs(new_r_z_min)<fabs(r[2]))*new_r_z_min  +  (fabs(new_r_z_plus)<fabs(new_r_z_min))*(fabs(new_r_z_plus)<fabs(r[2]))*new_r_z_plus    +   (fabs(r[2])<fabs(new_r_z_plus))*(fabs(r[2])<fabs(new_r_z_min))*r[2];

            distance = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
            index = distance/binSize;
            //cout<<index<<" "<<binSize<<" "<<numOfBins<<endl;
            if(index>maxIndex){
                maxIndex = index;
            }
            bins[index]+=1;

        }
    }
    cout<<"max: "<<maxIndex<<endl;
    return bins;
}

double CellSolver::gauss(double s, double mean){
    /*
      finds a random number from a gaussian distribution. the seed changes when the method is called
      */

    //return fmod(rand(),(4*s)) - 2*s;
    return r8_normal(mean,s,seed);
}
