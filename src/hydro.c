#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "chocobo_core.h"

CHRETURN get_total_energy(CHINT nel, CHFLOAT energy[nel], CHFLOAT rho[nel],
                          CHFLOAT u[nel], CHFLOAT v[nel], CHFLOAT mass[nel],
                          CHFLOAT **elwtc, CHINT **nodelist,
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
                                   CHINT **nodelist, CHFLOAT **ni,
                                   CHFLOAT **dndx, CHFLOAT **dndy,
                                   CHFLOAT **pdndx, CHFLOAT **pdndy,
                                   CHFLOAT **elwtc)
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
                         CHFLOAT **pdndx, CHFLOAT **pdndy,
                         CHINT **nodelist, CHFLOAT divvel[nel])
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


CHRETURN get_dt(CHINT nel, CHFLOAT rho[nel], CHFLOAT area[nel], CHFLOAT cc[nel],
                CHFLOAT q[nel], CHFLOAT time, CHFLOAT t0, CHFLOAT dtinit,
                CHFLOAT dtmax, CHFLOAT growth, CHFLOAT *dt, CHINT *dtcontrol)
{
    CHFLOAT dtmin = 1000.0;
    CHFLOAT dtold = *dt;
    CHFLOAT delta_t;

    *dtcontrol = 0;

    CELL_LOOP {
        delta_t = area[cell]/MAX(DENCUT,((pow(cc[cell],2))+2.0*(q[cell]/rho[cell])));
        delta_t = sqrt(delta_t)/2.0;

        if (area[cell] < 0.0) {
            printf("Negative area (%f) in cell %d\n", area[cell], cell+1);
        }
        if (delta_t < dtmin) {
            *dtcontrol = cell;
            dtmin = delta_t;
        }
    }

    if (time == t0) {
        *dt = MIN(dtmin, MIN(dtmax, dtinit));
    } else {
        *dt = MIN(dtmin, MIN(dtmax, growth*dtold));
    }

    if (time > 0.25 - *dt && time < 0.251005) {
        *dt = 0.25 - time;
    }
}


CHRETURN move_nodes(CHINT nnod, CHFLOAT dt, CHFLOAT x[nnod], CHFLOAT y[nnod],
                    CHFLOAT u[nnod], CHFLOAT v[nnod], CHFLOAT xout[nnod],
                    CHFLOAT yout[nnod])
{
    NODE_LOOP {
        xout[node] = x[node] + dt*u[node];
        yout[node] = y[node] + dt*v[node];
    }
}


CHRETURN calculate_volume(CHINT nel, CHINT nnod, CHFLOAT x[nnod],
                          CHFLOAT y[nnod], CHINT **nodelist,
                          CHFLOAT volume[nel], CHFLOAT area[nel])
{
    CHINT n1, n2, n3, n4;
    CHFLOAT a1, a2, a3;
    CHFLOAT b1, b2, b3;

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

        volume[cell] = 4.0*(a1*b3 - a3*b1);
        area[cell] = volume[cell];
    }
}


CHRETURN calculate_density(CHINT nel, CHFLOAT mass[nel], CHFLOAT volume[nel],
                           CHFLOAT rho[nel])
{
    CELL_LOOP {
        rho[cell] = mass[cell]/volume[cell];
    }
}


CHRETURN calculate_int_divv(CHINT zintdivvol, CHINT nel, CHINT nnod, CHFLOAT dt,
                            CHFLOAT vol[nel], CHFLOAT volold[nel],
                            CHFLOAT u[nnod], CHFLOAT v[nnod],
                            CHFLOAT **dndx, CHFLOAT **dndy,
                            CHINT **nodelist, CHFLOAT intdiv[nel])
{
    CHINT n;
    if (zintdivvol == 0) {
        CELL_LOOP {
            intdiv[cell] = 0.0;
        }
        CELL_LOOP {
            for (int j=0; j<4; j++) {
                n = nodelist[cell][j];
                intdiv[cell] += u[n]*dndx[cell][j] + v[n]*dndy[cell][j];
            }
        }
    } else {
        CELL_LOOP {
            intdiv[cell] = (vol[cell] - volold[cell])/dt;
        }
    }
}


CHRETURN calculate_energy(CHINT nel, CHFLOAT dt, CHFLOAT press[nel],
                          CHFLOAT visc[nel], CHFLOAT mass[nel],
                          CHFLOAT energy[nel], CHFLOAT intdiv[nel],
                          CHFLOAT enout[nel])
{
    CELL_LOOP {
        enout[cell] = energy[cell] - dt*(press[cell] + visc[cell])*intdiv[cell]/mass[cell];
    }
}


CHRETURN perfect_gas(CHINT nel, CHFLOAT energy[nel], CHFLOAT rho[nel],
                     CHINT material[nel], CHFLOAT gamma[MAXMAT],
                     CHFLOAT pressure[nel])
{
    CHINT mat;
    CELL_LOOP {
        mat = material[cell];
        pressure[cell] = (gamma[mat] - 1.0)*rho[cell]*energy[cell];
    }
}


CHRETURN momentum_calculation(CHINT nel, CHINT nnod, CHFLOAT dt, CHINT zantihg,
                              CHINT hgregtyp, CHFLOAT kappareg,
                              CHFLOAT u[nnod], CHFLOAT v[nnod], CHFLOAT x[nnod],
                              CHFLOAT y[nnod], CHFLOAT rho[nel],
                              CHFLOAT pressure[nel], CHFLOAT area[nel],
                              CHFLOAT cc[nel], CHFLOAT q[nel],
                              CHFLOAT **nint, CHFLOAT **dndx,
                              CHFLOAT **dndy, CHINT **nodelist,
                              CHINT znodbound[nel], CHFLOAT uout[nnod],
                              CHFLOAT vout[nnod])
{
    // Internal arrays
    CHFLOAT *massnod, *forcex, *forcey;

    // Other variables
    CHINT n;

    // Allocate these arrays
    massnod = ARRAY(CHFLOAT, nnod);
    forcex = ARRAY(CHFLOAT, nnod);
    forcey = ARRAY(CHFLOAT, nnod);

    // Set arrays to zero
    NODE_LOOP {
        massnod[node] = 0.0;
        forcex[node] = 0.0;
        forcey[node] = 0.0;
    }

    CELL_LOOP {
        for (int j=0; j<4; j++) {
            n = nodelist[cell][j];
            massnod[cell] += rho[cell]*nint[cell][j];
            forcex[cell] += (pressure[cell] + q[cell])*dndx[cell][j];
            forcey[cell] += (pressure[cell] + q[cell])*dndy[cell][j];
        }
    }

    NODE_LOOP {
        uout[node] = u[node] + dt*forcex[node]/massnod[node];
        vout[node] = v[node] + dt*forcey[node]/massnod[node];

        if (znodbound[node] == -1 || znodbound[node] == -3) uout[node] = u[node];
        if (znodbound[node] == -2 || znodbound[node] == -3) vout[node] = v[node];
    }

    // Free memory
    free(massnod);
    free(forcex);
    free(forcey);
}
