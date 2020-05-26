#ifndef CHOCOBO_MESHH
#define CHOCOBO_MESHH
CHRETURN get_mesh_size(char meshfile[CHSTRLEN], CHINT *nel, CHINT *nnod,
                       CHINT *nreg, CHINT *nmat);

CHRETURN read_mesh(char meshfile[CHSTRLEN], CHINT nel, CHINT nnod,
                   CHINT nodelist[nel][4], CHINT region[nel],
                   CHINT material[nel], CHINT regiontocell[nel][2],
                   CHFLOAT x[nnod], CHFLOAT y[nnod]);
#endif
