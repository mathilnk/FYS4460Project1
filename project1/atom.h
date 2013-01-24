#ifndef ATOM_H
#define ATOM_H

class Atom
{
public:
    Atom(double *position, double *velocity, char* element);
    double* position;
    double* velocity;
    char* element;

};

#endif // ATOM_H
