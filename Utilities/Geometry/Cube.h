#ifndef CUBE_H
#define CUBE_H

#include "Utilities/Geometry/Mesh.h"

class Cube : public Mesh {
public:
    Cube();
    Cube(const Cube& rhs);
    virtual ~Cube();

private:
    virtual void initFaces() override;
};

#endif // CUBE_H
