#ifndef GMRES_H
#define GMRES_H

#include "../CommonFiles/protos.h"
#include "../CommonFiles/Vector/Vector.h"

/* multiply a given CSR matrix to any Vector */
Vector matrix_vector_multiply_CSR(MAT* a, Vector b);

#endif
