#include "atom.h"

Atom::Atom(double* position, double* velocity, char* element)
{
    this->position = position;
    this->velocity = velocity;
    this->element = element;

}
