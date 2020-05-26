#include <stdio.h>
#include <stdlib.h>

#include "chocobo_core.h"

#define CHBUFFSIZE 255
#define ENTRYSIZEF 18

CHRETURN get_mesh_size(char meshfile[CHSTRLEN], CHINT *nel, CHINT *nnod,
                       CHINT *nreg, CHINT *nmat)
{
    FILE *fp;
    char buffer[CHBUFFSIZE];

    fp = fopen(meshfile, "r");

    fgets(buffer, CHBUFFSIZE, (FILE*)fp);
    *nel = atoi(buffer);

    fgets(buffer, CHBUFFSIZE, (FILE*)fp);
    *nnod = atoi(buffer);

    fgets(buffer, CHBUFFSIZE, (FILE*)fp);
    *nreg = atoi(buffer);
    *nmat = atoi(buffer);

    printf("Number of nodes: %d\n", *nnod);
    printf("Number of cells: %d\n", *nel);
    printf("Number of regions: %d\n", *nreg);
}

CHRETURN read_mesh(char meshfile[CHSTRLEN], CHINT nel, CHINT nnod,
                   CHINT nodelist[nel][4], CHINT region[nel],
                   CHINT material[nel], CHINT regiontocell[nel][2],
                   CHFLOAT x[nnod], CHFLOAT y[nnod])
{
    FILE *fp;
    char buffer[CHBUFFSIZE];

    fp = fopen(meshfile, "r");

    // Skip mesh sizes
    fgets(buffer, CHBUFFSIZE, (FILE*)fp);
    fgets(buffer, CHBUFFSIZE, (FILE*)fp);
    fgets(buffer, CHBUFFSIZE, (FILE*)fp);

    printf("%d\n", nel*ENTRYSIZEF);
}
