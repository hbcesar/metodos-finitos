/*
 *  Author: Josias Alexandre Oliveira
 *  emai: josiasalexandre@gmail.com
 *  copyleft =)
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Vector.h"

/* vector functions */
/* build a new vector with a given size */
Vector BuildVector(unsigned int s) {

    /* build a vector */
    Vector vector;

    /* assign the size value */
    vector.size = s;

    if (0 < s) {

        /* alloc the vector */
        vector.v = (double*) calloc(s, sizeof(double));
        if (NULL == vector.v) {

            printf("\nMemmory Error! Could not allocate the vector\n");
            exit(-20);

        }

    }

    /* return the vector */
    return(vector);

}

/* build a new vector with a given size */
Vector BuildVectorWithValue(unsigned int s, double value) {

    /* build a vector */
    Vector vector;

    /* helper */
    unsigned int i;

    /* assign the size value */
    vector.size = s;

    /* alloc the vector */
    vector.v = (double*) malloc(s*sizeof(double));
    if (NULL == vector.v) {

        printf("\nMemmory Error! Could not allocate the vector\n");
        exit(-21);

    }

    for (i = 0; i < s; i++) {

        vector.v[i] = value;

    }

    /* return the vector */
    return(vector);

}

/* build a new vector with random values */
Vector BuildVectorWithRandomValues(unsigned int s, unsigned int mod) {

    /* build a vector */
    Vector vector;

    /* helper */
    unsigned int i;

    /* assign the size value */
    vector.size = s;

    /* alloc the vector */
    vector.v = (double*) malloc(s*sizeof(double));
    if (NULL == vector.v) {

        printf("\nMemmory Error! Could not allocate the vector\n");
        exit(-21);

    }

    for (i = 0; i < s; i++) {

        vector.v[i] = rand() % mod;

    }

    /* return the vector */
    return(vector);
}

/* delete a given vector */
void DeleteVector(Vector vector) {

    /* remove the vector */
    free(vector.v);

    return;

}

/* copy a vector */
Vector CopyVector(Vector vector) {

    /* helpers */
    unsigned int i;
    Vector newVector;

    /* build a new vector with the same size */
    newVector = BuildVector(vector.size);

    /* copy each element */
    for (i = 0; i < vector.size; i++) {

        /* copy the element */
        newVector.v[i] = vector.v[i];

    }

    /* return the new vector */
    return(newVector);

}

/* Copy the vector A to the vector B */
void CopyVectorAToB(Vector *a, Vector *b) {

    unsigned int i;

    for (i = 0; i < a->size; i++) {

        /* copy the current element */
        b->v[i] = a->v[i];

    }

    return;

}

/* inner product */
double InnerProduct(Vector a, Vector b) {

    /* helpers */
    unsigned int i, size = a.size;
    double result = 0.0f;

    if (size != b.size) {

        printf("Error! The vectors does not have the same sizes!\n");
        exit(-22);

    } else {

        /* multiply each element */
        for (i = 0; i < size; i++) {

            /* update the result */
            result += a.v[i]*b.v[i];

        }

    }

    /* return the desired result */
    return(result);

}

/* cross product */
Vector CrossProduct(Vector a, Vector b) {

    /* build the new vector */
    Vector newVector = BuildVector(3);

    if (3 != a.size || 3 != b.size) {

        /* error! */
        printf("Error! The input vector are not 3D vectors!\n");

    } else {

        /* get the cross product */
        newVector.v[0] = a.v[1]*b.v[2] - a.v[2]*b.v[1];
        newVector.v[1] = a.v[2]*b.v[0] - a.v[0]*b.v[2];
        newVector.v[2] = a.v[0]*b.v[1] - a.v[1]*b.v[0];

    }

    /* return the new Vector */
    return newVector;

}

/* Euclidean norm */
double EuclideanNorm(Vector vector) {

    /* get the vector size */
    unsigned int v_size = vector.size, i;

    /* the resulting norm */
    double norm = 0;

    if (0 < v_size) {

        /* syntatic sugvectorr */
        double *v = vector.v;

        for (i = 0; i < v_size; i++) {

            norm += v[i]*v[i];

        }

    } else {

        return 0.0;
    }

    return sqrt(norm);

}

/* infinity norm */
double InfinityNorm(Vector vector) {

    /* get the vector size */
    unsigned int v_size = vector.size, i;

    /* syntatic sugar */
    double *v = vector.v;

    /* the resulting norm */
    double norm = fabs(v[0]);

    /* the potential new norm */
    double p;

    for (i = 1; i < v_size; i++) {

        p = fabs(v[i]);

        if (norm < p) {

            norm = p;

        }

    }

    return norm;

}

/* normalization */
void ScaleVector(Vector vector, double value) {

    /* get the vector size */
    unsigned int v_size = vector.size, i;

    /* syntatic sugvectorr */
    double *v = vector.v;

    for (i = 0; i < v_size; i++) {

        v[i] *= value;

    }

}

/* show a vector */
void ShowVector(Vector vector) {

    unsigned int i;

    printf("[ ");
    for (i = 0; i < vector.size; i++) {

        printf("%.5f ", vector.v[i]);

    }
    printf("]\n");

}
