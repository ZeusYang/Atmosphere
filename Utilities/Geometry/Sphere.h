#ifndef SPHERE_H
#define SPHERE_H

#include "Utilities/Geometry/Mesh.h"

class Sphere : public Mesh {
private:
    int longitude,latitude;//经度间隔，维度间隔
    virtual void initFaces();

public:
    Sphere(int m,int n);
    Sphere(const Sphere& rhs);
    virtual ~Sphere();
};

#endif // SPHERE_H
