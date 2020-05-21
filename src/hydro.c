#include <math.h>

#include "chocobo_core.h"

CHRETURN get_total_energy(CHINT nel, CHFLOAT energy[nel], CHFLOAT rho[nel],
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

    CELL_LOOP {
        tek = 0.0;
        VERTEX_LOOP {
            node = nodelist[cell][vertex];
            tek += 0.5*rho[cell]*elwtc[cell][vertex]*(
                pow(u[node],2) + pow(v[node],2)
            )*TWOPI;
        }
        element_energy = mass[cell]*energy[cell]*TWOPI + tek;
        *total_ke += tek;
        *total_ie += mass[cell]*energy[cell]*TWOPI;
        *total_energy += element_energy;
    }
}


CHRETURN calculate_mass(
    CHINT nel, CHFLOAT volume[nel], CHFLOAT rho[nel], CHFLOAT mass[nel]
)
{
    CELL_LOOP {
        mass[cell] = volume[cell]*rho[cell];
    }
}


CHRETURN calculate_finite_elements(CHINT nel, CHFLOAT x[nel], CHFLOAT y[nel],
                                   CHINT nodelist[nel][4], CHFLOAT ni[nel][4],
                                   CHFLOAT dndx[nel][4], CHFLOAT dndy[nel][4],
                                   CHFLOAT pdndx[nel][4], CHFLOAT pdndy[nel][4],
                                   CHFLOAT elwtc[nel][4])
{
    CHFLOAT a1, a2, a3;
    CHFLOAT b1, b2, b3;
    CHFLOAT jacobian;
    CHINT n1, n2, n3, n4;

    CELL_LOOP {
        n1 = nodelist[cell][0];
        n2 = nodelist[cell][1];
        n3 = nodelist[cell][2];
        n4 = nodelist[cell][3];

        a1 = 0.25*(-x[n1] + x[n2] + x[n3] - x[n4]);
        a2 = 0.25*(x[n1] - x[n2] + x[n3] - x[n4]);
        a3 = 0.25*(-x[n1] - x[n2] + x[n3] + x[n4]);

        b1 = 0.25*(-y[n1] + y[n2] + y[n3] - y[n4]);
        b2 = 0.25*(y[n1] - y[n2] + y[n3] - y[n4]);
        b3 = 0.25*(-y[n1] - y[n2] + y[n3] + y[n4]);

        // Used for momentum update
        ni[cell][0] = ((3.0*b3 - b2)*(3.0*a1 - a2) - (3.0*a3 - a2)*(3.0*b1 - b2))/9.0;
        ni[cell][1] = ((3.0*b3 + b2)*(3.0*a1 - a2) - (3.0*a3 + a2)*(3.0*b1 - b2))/9.0;
        ni[cell][2] = ((3.0*b3 + b2)*(3.0*a1 + a2) - (3.0*a3 + a2)*(3.0*b1 + b2))/9.0;
        ni[cell][3] = ((3.0*b3 - b2)*(3.0*a1 + a2) - (3.0*a3 - a2)*(3.0*b1 + b2))/9.0;

        // Integrated d(dndx) and d(dndy) for intdiv and energy
        dndx[cell][0] = -b3 + b1;
        dndx[cell][1] = b3 + b1;
        dndx[cell][2] = b3 - b1;
        dndx[cell][3] = -b3 - b1;

        dndy[cell][0] = a3 - a1;
        dndy[cell][1] = -a3 - a1;
        dndy[cell][2] = -a3 + a1;
        dndy[cell][3] = a3 + a1;

        // Jacobian and d(dndx), d(dndy) for div and viscosity
        jacobian = a1*b3 - a3*b1;
        VERTEX_LOOP {
            pdndx[cell][vertex] = 0.25*dndx[cell][vertex]/jacobian;
            pdndy[cell][vertex] = 0.25*dndy[cell][vertex]/jacobian;
        }

        // Copy ni into elwtc
        VERTEX_LOOP {
            elwtc[cell][vertex] = ni[cell][vertex];
        }
    }
}


CHRETURN calculate_div_v(CHINT nel, CHFLOAT u[nel], CHFLOAT v[nel],
                         CHFLOAT pdndx[nel][4], CHFLOAT pdndy[nel][4],
                         CHINT nodelist[nel][4], CHFLOAT divvel[nel])
{
    CHINT n;

    CELL_LOOP {
        divvel[cell] = 0.0;
    }

    CELL_LOOP {
        VERTEX_LOOP {
            n = nodelist[cell][vertex];
            divvel[cell] += u[n]*pdndx[cell][vertex] + v[n]*pdndy[cell][vertex];
        }
    }
}


CHRETURN calculate_soundspeed(CHINT nel, CHFLOAT pressure[nel], CHFLOAT rho[nel],
                              CHINT material[nel], CHFLOAT gamma[MAXMAT],
                              CHFLOAT soundspeed[nel])
{
    CHINT mat;

    CELL_LOOP {
        mat = material[cell];
        soundspeed[cell] = sqrt(gamma[mat]*pressure[cell]/rho[cell]);
    }
}

CHRETURN calculate_q(CHINT nel, CHFLOAT rho[nel], CHFLOAT soundspeed[nel],
                     CHFLOAT divvel[nel], CHFLOAT area[nel], CHFLOAT cq,
                     CHFLOAT cl, CHFLOAT q[nel])
{
    CHFLOAT dudx;

    CELL_LOOP {
        if (divvel[cell] < 0.0) {
            dudx = sqrt(area[cell])*divvel[cell];
            q[cell] = cq*rho[cell]*dudx*dudx +
                      cl*rho[cell]*soundspeed[cell]*fabs(dudx);
        } else {
            q[cell] = 0.0;
        }
    }
}
