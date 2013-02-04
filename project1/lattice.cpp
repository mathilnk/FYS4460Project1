#include "lattice.h"




Lattice::Lattice(){}
Lattice::Lattice(int Nx, int Ny, int Nz, string element, double b, double T, double mass)
{
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;
    this->numberOfAtoms = N*4;
    this->element = element;
    this->b = b;
    this->T = T;
    this->k =1.38e-23; //J/K
    this->mass = mass;
    this->s = sqrt(k*T/mass);
    //this->s = 144.0;
    cout<<s<<endl;
    this->mean = 0;
    this->pi = 3.1415;
    this->seed = 123456789;


    makeLattice();
}

void Lattice::makeLattice(){
    vec posBase(3);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz; k++){
                posBase<<i*b<<j*b<<k*b;
                findPosAndMakeAtoms2(posBase);
            }
        }
    }
    //makeEndAtoms();
}


void Lattice::makeEndAtoms(){
//top square
    vec pos, vec;
    Atom *atom1,*atom2,*atom3,*atom4,*atom5,*atom6,*atom7;
    vec<<0<<0<<0;
    pos<<(Nx-1)*b<<(Ny-1)*b<<Nz*b;
    atom1= new Atom(pos,vec, element);
    allAtoms.push_back(atom1);

    pos[0] +=b;
    atom2= new Atom(pos,vec, element);

    allAtoms.push_back(atom2);

    pos[1]+=b;
    atom3= new Atom(pos,vec,element);

    allAtoms.push_back(atom3);


    pos[0]-=b;
    atom4= new Atom(pos,vec, element);

    allAtoms.push_back(atom4);

    //bottom square

    pos<<Nx*b<<(Ny-1)*b<<(Nz-1)*b;
    atom5= new Atom(pos,vec, element);

    allAtoms.push_back(atom5);

    pos[1] +=b;
    atom6= new Atom(pos,vec, element);

    allAtoms.push_back(atom6);

    pos[0]-=b;
    atom7= new Atom(pos,vec, element);

    allAtoms.push_back(atom7);







}

void Lattice::writeVMDfile(const char *Filename, string comment){
    ofstream writeToFile;
    numberOfAtoms = allAtoms.size();
    writeToFile.open(Filename);
    writeToFile<<numberOfAtoms<<endl;
    writeToFile<<comment<<endl;
    for(int i=0;i<numberOfAtoms;i++){
        Atom* a = allAtoms[i];
        writeToFile<<element<<" ";
        for(int j=0;j<3;j++){
            writeToFile<<a->position[j]<<" ";
        }
        for(int k=0;k<3;k++){
            writeToFile<<a->velocity[k]<<" ";
        }
        writeToFile<<"\n";
    }


}

void Lattice::findPosAndMakeAtoms2(vec posBase){
    ////////////////////
    //modifiser denne//
    //////////////////
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
    allAtoms.push_back(xy);
    allAtoms.push_back(yz);
    allAtoms.push_back(zx);
    allAtoms.push_back(base);

}
void Lattice::findPosAndMakeAtoms(vec posBase){
    vec posxy(3);
    vec posyz(3);
    vec poszx(3);
    posxy<<0.5*b<<0.5*b<<0;
    posyz<<0<<0.5*b<<0.5*b;
    poszx<<0.5*b<<0<<0.5*b;
    posxy+=posBase;
    posyz+=posBase;
    poszx+=posBase;
    vec velo(3);
    velo<<0<<0<<0;
    Atom *xy = new Atom(posxy,velo, element);
    Atom *yz = new Atom(posyz,velo, element);
    Atom *zx = new Atom(poszx,velo, element);
    Atom *base = new Atom(posBase,velo, element);
    allAtoms.push_back(xy);
    allAtoms.push_back(yz);
    allAtoms.push_back(zx);
    allAtoms.push_back(base);
}

double Lattice::gauss(double s, double mean){

    return r8_normal(mean,s,seed);
}



