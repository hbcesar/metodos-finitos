/*
 *  Author: Josias Alexandre Oliveira
 *  emai: josiasalexandre@gmail.com
 *  copyleft =)
 */

#ifndef VECTOR_DATATYPE_H
#define VECTOR_DATATYPE_H

typedef struct Vector {

    /* how many elements? */
    unsigned int size;

    /* the actual vector */
    double *v;

} Vector;

/* vector functions */
/* build a new vector with a given size */
Vector BuildVector(unsigned int s);

/* build a new vector with a given size and value */
Vector BuildVectorWithValue(unsigned int s, double value);

/* build a new vector with random values */
Vector BuildVectorWithRandomValues(unsigned int s, unsigned int mod);

/* delete a given vector */
void DeleteVector(Vector vector);

/* copy a vector */
Vector CopyVector(Vector vector);

/* Copy the vector A to the vector B */
void CopyVectorAToB(Vector *a, Vector *b);

/* inner product */
double InnerProduct(Vector a, Vector b);

/* cross product */
Vector CrossProduct(Vector a, Vector b);

/* Euclidean norm */
double EuclideanNorm(Vector a);

/* infinity norm */
double InfinityNorm(Vector vector);

/* normalization */
inline void ScaleVector(Vector vector, double value);

/* show a vector */
void ShowVector(Vector vector);

#endif
