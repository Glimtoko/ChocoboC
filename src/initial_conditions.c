#include "chocobo_core.h"

#include <stdio.h>

CHRETURN set_initial_conditions(char eosfile[CHSTRLEN], CHINT nel, CHINT nmat,
                                CHINT material[nel], CHFLOAT pre[nel],
                                CHFLOAT rho[nel], CHFLOAT gamma[MAXMAT])
{
    FILE *fp = fopen(eosfile, "r");

    CHFLOAT data[3];

    MATERIAL_LOOP {
        for (int i=0; i<3; i++) {
            fscanf(fp, "%lf", &data[i]);
        }
        printf("\nMaterial %d\n", mat+1);
        printf("Initial density: %f\n", data[0]);
        printf("Initial pressure: %f\n", data[1]);
        printf("Gamma: %f\n", data[2]);

        CELL_LOOP {
            if (material[cell] == mat+1) {
                rho[cell] = data[0];
                pre[cell] = data[1];
            }
        }
        gamma[mat] = data[2];
    }
}
