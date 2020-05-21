#include <math.h>

#include "chocobo_core.h"

void get_total_energy(CHINT nel, CHFLOAT energy[nel], CHFLOAT rho[nel],
                      CHFLOAT u[nel], CHFLOAT v[nel], CHFLOAT mass[nel],
                      CHFLOAT elwtc[nel][4], CHINT nodelist[nel][4],
                      CHFLOAT *total_energy, CHFLOAT *total_ke,
                      CHFLOAT *total_ie)
{
    *total_energy = 0.0;
    *total_ie = 0.0;
    *total_ke = 0.0;

    CHFLOAT tek;
    CHFLOAT element_energy;
    CHINT node;

    for (int i=0; i<nel; i++) {
        tek = 0.0;
        for (int j=0; j<4; j++) {
            node = nodelist[i][j];
            tek += 0.5*rho[i]*elwtc[i][j]*(pow(u[node],2) + pow(v[node],2))*TWOPI;
        }
        element_energy = mass[i]*energy[i]*TWOPI + tek;
        *total_ke += tek;
        *total_ie += mass[i]*energy[i]*TWOPI;
        *total_energy += element_energy;
    }
}


void calculate_mass(
    CHINT nel, CHFLOAT volume[nel], CHFLOAT rho[nel], CHFLOAT mass[nel]
)
{
    for (CHINT i=0; i<nel; i++) {
        mass[i] = volume[i]*rho[i];
    }
}


void calculate_finite_elements(CHINT nel, CHFLOAT x[nel], CHFLOAT y[nel],
                               CHINT nodelist[nel][4], CHFLOAT ni[nel][4],
                               CHFLOAT dndx[nel][4], CHFLOAT dndy[nel][4],
                               CHFLOAT pdndx[nel][4], CHFLOAT pdndy[nel][4],
                               CHFLOAT elwtc[nel][4])
{
    CHFLOAT a1, a2, a3;
    CHFLOAT b1, b2, b3;
    CHFLOAT jacobian;
    CHINT n1, n2, n3, n4;

    for (int i=0; i<nel; i++) {
        n1 = nodelist[i][0];
        n2 = nodelist[i][1];
        n3 = nodelist[i][2];
        n4 = nodelist[i][3];

        a1 = 0.25*(-x[n1] + x[n2] + x[n3] - x[n4]);
        a2 = 0.25*(x[n1] - x[n2] + x[n3] - x[n4]);
        a3 = 0.25*(-x[n1] - x[n2] + x[n3] + x[n4]);

        b1 = 0.25*(-y[n1] + y[n2] + y[n3] - y[n4]);
        b2 = 0.25*(y[n1] - y[n2] + y[n3] - y[n4]);
        b3 = 0.25*(-y[n1] - y[n2] + y[n3] + y[n4]);

        // Used for momentum update
        ni[i][0] = ((3.0*b3 - b2)*(3.0*a1 - a2) - (3.0*a3 - a2)*(3.0*b1 - b2))/9.0;
        ni[i][1] = ((3.0*b3 + b2)*(3.0*a1 - a2) - (3.0*a3 + a2)*(3.0*b1 - b2))/9.0;
        ni[i][2] = ((3.0*b3 + b2)*(3.0*a1 + a2) - (3.0*a3 + a2)*(3.0*b1 + b2))/9.0;
        ni[i][3] = ((3.0*b3 - b2)*(3.0*a1 + a2) - (3.0*a3 - a2)*(3.0*b1 + b2))/9.0;

        // Integrated d(dndx) and d(dndy) for intdiv and energy
        dndx[i][0] = -b3 + b1;
        dndx[i][1] = b3 + b1;
        dndx[i][2] = b3 - b1;
        dndx[i][3] = -b3 - b1;

        dndy[i][0] = a3 - a1;
        dndy[i][1] = -a3 - a1;
        dndy[i][2] = -a3 + a1;
        dndy[i][3] = a3 + a1;

        // Jacobian and d(dndx), d(dndy) for div and viscosity
        jacobian = a1*b3 - a3*b1;
        for (int j=0; j<4; j++) {
            pdndx[i][j] = 0.25*dndx[i][j]/jacobian;
            pdndy[i][j] = 0.25*dndy[i][j]/jacobian;
        }

        // Copy ni into elwtc
        for (int j=0; j<4; j++) {
            elwtc[i][j] = ni[i][j];
        }
    }
}
