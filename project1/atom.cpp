#include "atom.h"

Atom::Atom(){
}

Atom::Atom(vec position, vec velocity, string element)
{
    /*
      Constructor
      */
    this->position = position;
    this->velocity = velocity;
    this->element = element;

}
