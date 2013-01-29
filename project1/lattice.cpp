#include "lattice.h"





Lattice::Lattice(int Nx, int Ny, int Nz, string element, double b)
{
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;
    this->numberOfAtoms = N*4;
    this->element = element;
    this->b = b;
    makeLattice();
}

void Lattice::makeLattice(){
    vec posBase(3);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz; k++){
                posBase<<i*b<<j*b<<k*b;
                findPosAndMakeAtoms(posBase);
            }
        }
    }
}



void Lattice::writeVMDfile(const char *Filename, string comment){
    ofstream writeToFile;
    writeToFile.open(Filename);
    writeToFile<<numberOfAtoms<<endl;
    writeToFile<<comment<<endl;
    for(int i=0;i<numberOfAtoms;i++){
        Atom* a = allAtoms[i];
        writeToFile<<element;
        for(int j=0;j<6;j++){
            writeToFile<<a->position[j];
        }
        writeToFile<<"\n";
    }


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
