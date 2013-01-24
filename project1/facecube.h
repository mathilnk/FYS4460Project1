#ifndef FACECUBE_H
#define FACECUBE_H

class FaceCube
{
public:
    FaceCube(double b, char *element);
    double b;
    char* element;
    double* origo, *xy, *yz, *zx;
};

#endif // FACECUBE_H
