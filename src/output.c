#include "chocobo_core.h"

#include <string.h>
#include <stdio.h>

CHRETURN output(char baseloc[CHSTRLEN], CHINT step, CHINT nel, CHINT nnod,
                CHFLOAT time, CHFLOAT xv[nnod], CHFLOAT yv[nnod],
                CHFLOAT uv[nnod], CHFLOAT vv[nnod], CHINT **nodelist,
                CHFLOAT rho[nel], CHFLOAT pre[nel], CHFLOAT en[nel],
                CHFLOAT volume[nel], CHINT *filenum)
{
    char filename[CHSTRLEN];
    FILE *fp, *fp2;

    // Nodelist
    sprintf(filename, "%s/results/nodelist_no%03d.txt", baseloc, *filenum);
    fp = fopen(filename, "w");
    CELL_LOOP {
        fprintf(fp, "%d %d %d %d\n", nodelist[cell][0]+1, nodelist[cell][1]+1,
                nodelist[cell][2]+1, nodelist[cell][3]+1);
    }
    fclose(fp);

    // Cell-centred quantities
    sprintf(filename, "%s/results/mypre_no%03d.txt", baseloc, *filenum);
    fp = fopen(filename, "w");
    sprintf(filename, "%s/results/vol_no%03d.txt", baseloc, *filenum);
    fp2 = fopen(filename, "w");
    CELL_LOOP {
        fprintf(fp, "%d %f %f %f\n", cell, pre[cell], rho[cell], en[cell]);
        fprintf(fp2, "%d %f\n", cell, volume[cell]);
    }
    fclose(fp);
    fclose(fp2);

    // Node-centred quantities
    sprintf(filename, "%s/results/myvel_no%03d.txt", baseloc, *filenum);
    fp = fopen(filename, "w");
    sprintf(filename, "%s/results/myx_no%03d.txt", baseloc, *filenum);
    fp2 = fopen(filename, "w");
    NODE_LOOP {
        fprintf(fp, "%f %f %f\n", time, uv[node], vv[node]);
        fprintf(fp2, "%f %f\n", xv[node], yv[node]);
    }
    fclose(fp);
    fclose(fp2);

    filenum++;
}
