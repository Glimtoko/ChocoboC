#ifndef CHOCOBOC_DATAH
#define CHOCOBOC_DATAH

#include "chocobo_core.h"

// Main struct
struct MeshT {
    CHINT nel;
    CHINT nnod;
    CHINT nreg;
    CHINT nmat;

    // "logicals" - actually integers
    CHINT *znodbound;

    //       ..arrays..
    CHINT **nodelist;
    CHINT *region;
    CHINT *material;
    CHINT *nodetype;
    CHINT **regiontocell;

    CHFLOAT *gamma;

    CHFLOAT *xv;
    CHFLOAT *yv;
    CHFLOAT *pre;
    CHFLOAT *rho;
    CHFLOAT *en;
    CHFLOAT *cc;
    CHFLOAT *qq;
    CHFLOAT *massel;
    CHFLOAT *volel;
    CHFLOAT *volelold;
    CHFLOAT *area;
    CHFLOAT *uv;
    CHFLOAT *vv;

    CHFLOAT *uvold;
    CHFLOAT *vvold;
    CHFLOAT *uvbar;
    CHFLOAT *vvbar;
    CHFLOAT *xv05;
    CHFLOAT *yv05;
    CHFLOAT *pre05;
    CHFLOAT *rho05;
    CHFLOAT *en05;
    CHFLOAT *volel05;
    CHFLOAT *divint;
    CHFLOAT *divvel;

    // FEM
    CHFLOAT **nint;
    CHFLOAT **elwtc;
    CHFLOAT **dndx;
    CHFLOAT **dndy;
    CHFLOAT **pdndx;
    CHFLOAT **pdndy;
};

struct InputT {
    CHFLOAT cq;
    CHFLOAT cl;

    CHFLOAT dtmax;
    CHFLOAT dtinit;
    CHFLOAT growth;

    CHFLOAT t0;
    CHFLOAT tend;
    CHFLOAT tout;

    CHINT zintdivvol;

    CHINT stepcount;
};

#endif
