#ifndef CHOCOBO_MESHH
#define CHOCOBO_MESHH
CHRETURN get_mesh_size(char meshfile[CHSTRLEN], CHINT *nel, CHINT *nnod,
                       CHINT *nreg, CHINT *nmat);

CHRETURN read_mesh(char meshfile[CHSTRLEN], CHINT nel, CHINT nnod,
                   CHINT **nodelist, CHINT nodetype[nel], CHINT region[nel],
                   CHINT material[nel], CHINT **regiontocell,
                   CHFLOAT x[nnod], CHFLOAT y[nnod]);
#endif
