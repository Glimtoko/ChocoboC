#include <stdio.h>
#include <stdlib.h>

#include "chocobo_core.h"

#define CHBUFFSIZE 255
#define ENTRYSIZEF 18

CHRETURN get_mesh_size(char meshfile[CHSTRLEN], CHINT *nel, CHINT *nnod,
                       CHINT *nreg, CHINT *nmat)
{
    FILE *fp = fopen(meshfile, "r");

    fscanf(fp, "%d", nnod);
    fscanf(fp, "%d", nel);
    fscanf(fp, "%d", nreg);
    *nmat = *nreg;

    printf("Number of nodes: %d\n", *nnod);
    printf("Number of cells: %d\n", *nel);
    printf("Number of regions: %d\n", *nreg);
}

CHRETURN read_mesh(char meshfile[CHSTRLEN], CHINT nel, CHINT nnod,
                   CHINT **nodelist, CHINT nodetype[nel], CHINT region[nel],
                   CHINT material[nel], CHINT **regiontocell,
                   CHFLOAT x[nnod], CHFLOAT y[nnod])
{
    FILE *fp = fopen(meshfile, "r");

    // Skip mesh sizes
    CHINT dummy;
    CHINT nreg;

    fscanf(fp, "%d", &dummy);
    fscanf(fp, "%d", &dummy);
    fscanf(fp, "%d", &nreg);


    // X coordinates
    NODE_LOOP {
        fscanf(fp, "%le", &x[node]);
    }

    // Y coordinates
    NODE_LOOP {
        fscanf(fp, "%le", &y[node]);
    }

    // Node type
    NODE_LOOP {
        fscanf(fp, "%d", &nodetype[node]);
    }

    // Region list
    CELL_LOOP {
        fscanf(fp, "%d", &region[cell]);
    }

    // Material list
    CELL_LOOP {
        fscanf(fp, "%d", &material[cell]);
    }

    // Node list - requires reordering
    CHINT *nodelistraw = malloc(nel * 4 * sizeof(int));
    for (int i=0; i<nel*4; i++) {
        fscanf(fp, "%d", &nodelistraw[i]);
    }

    // Store in nodelist
    int j = 0;
    CELL_LOOP {
        VERTEX_LOOP {
            nodelist[cell][vertex] = nodelistraw[j] - 1;
            j++;
        }
    }

    // Free the temporary array
    free(nodelistraw);

    // Region to cell iterators
    CHINT d[3];
    REGION_LOOP {
        for (int i=0; i<3; i++) {
            fscanf(fp, "%d", &d[i]);
            regiontocell[reg][0] = d[1];
            regiontocell[reg][1] = d[2];
        }
    }

}
