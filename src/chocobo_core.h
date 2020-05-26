#ifndef CHOCOBOC_COREH
#define CHOCOBOC_COREH

// Data types
#define CHINT int
#define CHFLOAT double
#define CHRETURN void
#define CHSTRLEN 100

// Array limits
#define MAXMAT 10

// Pi
#define PI 3.14159265358979324
#define TWOPI 2.0*3.14159265358979324

// Cutoffs
#define DTMINHG 0.000000008
#define DENCUT 1.000000000E-6

// LOOPS
#define CELL_LOOP for (int cell=0; cell<nel; cell++)
#define NODE_LOOP for (int node=0; node<nnod; node++)
#define VERTEX_LOOP for (int vertex=0; vertex<4; vertex++)

// Useful functions
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// Allocation "functions"
#define ARRAY(a,n)  (a*) malloc(n * sizeof(a))
// #define SETARRAY2D(t,r,c)  (t*) malloc(r * sizeof(t)); for(int ii=0;i<r;i++){}
#define SETARRAY2D(array,t,r,c) \
    array = malloc(r*sizeof(t*)); \
    for (int ii=0;ii<r;ii++) { \
        array[ii] = malloc(c*sizeof(int)); \
    }

#endif
