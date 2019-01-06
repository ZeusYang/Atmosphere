#ifndef QUAD_H
#define QUAD_H

#include "Utilities/Geometry/Mesh.h"

class Quad : public Mesh {
private:
    virtual void initFaces() override;
public:
    Quad();
    Quad(const Quad& rhs);
    virtual ~Quad();
};

#endif // QUAD_H
