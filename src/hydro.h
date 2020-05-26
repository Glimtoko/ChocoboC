#ifndef CHOCOBOC_HYDROH
#define CHOCOBOC_HYDROH

#include "chocobo_core.h"

CHRETURN get_total_energy(CHINT nel, CHFLOAT energy[nel], CHFLOAT rho[nel],
                          CHFLOAT u[nel], CHFLOAT v[nel], CHFLOAT mass[nel],
                          CHFLOAT elwtc[nel][4], CHINT nodelist[nel][4],
                          CHFLOAT *total_energy, CHFLOAT *total_ke,
                          CHFLOAT *total_ie);


CHRETURN calculate_mass(CHINT nel, CHFLOAT volume[nel], CHFLOAT rho[nel],
                        CHFLOAT mass[nel]);


CHRETURN calculate_finite_elements(CHINT nel, CHFLOAT x[nel], CHFLOAT y[nel],
                                   CHINT nodelist[nel][4], CHFLOAT ni[nel][4],
                                   CHFLOAT dndx[nel][4], CHFLOAT dndy[nel][4],
                                   CHFLOAT pdndx[nel][4], CHFLOAT pdndy[nel][4],
                                   CHFLOAT elwtc[nel][4]);


CHRETURN calculate_div_v(CHINT nel, CHFLOAT u[nel], CHFLOAT v[nel],
                         CHFLOAT pdndx[nel][4], CHFLOAT pdndy[nel][4],
                         CHINT nodelist[nel][4], CHFLOAT divvel[nel]);


CHRETURN calculate_soundspeed(CHINT nel, CHFLOAT pressure[nel], CHFLOAT rho[nel],
                              CHINT material[nel], CHFLOAT gamma[MAXMAT],
                              CHFLOAT soundspeed[nel]);


CHRETURN calculate_q(CHINT nel, CHFLOAT rho[nel], CHFLOAT soundspeed[nel],
                     CHFLOAT divvel[nel], CHFLOAT area[nel], CHFLOAT cq,
                     CHFLOAT cl, CHFLOAT q[nel]);


CHRETURN get_dt(CHINT nel, CHFLOAT rho[nel], CHFLOAT area[nel], CHFLOAT cc[nel],
                CHFLOAT q[nel], CHFLOAT time, CHFLOAT t0, CHFLOAT dtinit,
                CHFLOAT dtmax, CHFLOAT growth, CHFLOAT *dt, CHINT *dtcontrol);


CHRETURN move_nodes(CHINT nnod, CHFLOAT dt, CHFLOAT x[nnod], CHFLOAT y[nnod],
                    CHFLOAT u[nnod], CHFLOAT v[nnod], CHFLOAT xout[nnod],
                    CHFLOAT yout[nnod]);


CHRETURN calculate_density(CHINT nel, CHFLOAT mass[nel], CHFLOAT volume[nel],
                           CHFLOAT rho[nel]);


CHRETURN calculate_int_divv(CHINT zintdivvol, CHINT nel, CHINT nnod, CHFLOAT dt,
                            CHFLOAT vol[nel], CHFLOAT volold[nel],
                            CHFLOAT u[nnod], CHFLOAT v[nnod],
                            CHFLOAT dndx[nel][4], CHFLOAT dndy[nel][4],
                            CHINT nodelist[nel][4], CHFLOAT intdiv[nel]);


CHRETURN calculate_energy(CHINT nel, CHFLOAT dt, CHFLOAT press[nel],
                          CHFLOAT visc[nel], CHFLOAT mass[nel],
                          CHFLOAT energy[nel], CHFLOAT intdiv[nel],
                          CHFLOAT enout[nel]);


CHRETURN perfect_gas(CHINT nel, CHFLOAT energy[nel], CHFLOAT rho[nel],
                     CHINT material[nel], CHFLOAT gamma[MAXMAT],
                     CHFLOAT pressure[nel]);


CHRETURN momentum_calculation(CHINT nel, CHINT nnod, CHFLOAT dt, CHINT zantihg,
                              CHINT hgregtyp, CHFLOAT kappareg,
                              CHFLOAT u[nnod], CHFLOAT v[nnod], CHFLOAT x[nnod],
                              CHFLOAT y[nnod], CHFLOAT rho[nel],
                              CHFLOAT pressure[nel], CHFLOAT area[nel],
                              CHFLOAT cc[nel], CHFLOAT q[nel],
                              CHFLOAT nint[nel][4], CHFLOAT dndx[nel][4],
                              CHFLOAT dndy[nel][4], CHINT nodelist[nel][4],
                              CHINT znodbound[nel], CHFLOAT uout[nnod],
                              CHFLOAT vout[nnod]);

#endif
