#include "chocobo_core.h"
#include "hydro.h"
#include "data.h"
#include "mesh.h"
#include "initial_conditions.h"

#include <stdlib.h>
#include <stdio.h>

int set_store(struct MeshT *mesh, CHINT nel, CHINT nnod, CHINT nreg,
             CHINT nmat)
{
    mesh->nel = nel;
    mesh->nnod = nnod;
    mesh->nreg = nreg;
    mesh->nmat = nmat;

    mesh->znodbound = ARRAY(CHINT, nnod);

    SETARRAY2D(mesh->nodelist, CHINT, nel, 4);
    mesh->region = ARRAY(CHINT, nel);
    mesh->material = ARRAY(CHINT, nel);
    mesh->nodetype = ARRAY(CHINT, nnod);
    SETARRAY2D(mesh->regiontocell, CHINT, nel, 2);

    mesh->gamma = ARRAY(CHFLOAT, MAXMAT);

    mesh->xv = ARRAY(CHFLOAT, nnod);
    mesh->yv = ARRAY(CHFLOAT, nnod);
    mesh->pre = ARRAY(CHFLOAT, nel);
    mesh->rho = ARRAY(CHFLOAT, nel);
    mesh->en = ARRAY(CHFLOAT, nel);
    mesh->cc = ARRAY(CHFLOAT, nel);
    mesh->qq = ARRAY(CHFLOAT, nel);
    mesh->massel = ARRAY(CHFLOAT, nel);
    mesh->volel = ARRAY(CHFLOAT, nel);
    mesh->volelold = ARRAY(CHFLOAT, nel);
    mesh->area = ARRAY(CHFLOAT, nel);
    mesh->uv = ARRAY(CHFLOAT, nnod);
    mesh->vv = ARRAY(CHFLOAT, nnod);

    mesh->uvold = ARRAY(CHFLOAT, nnod);
    mesh->vvold = ARRAY(CHFLOAT, nnod);
    mesh->uvbar = ARRAY(CHFLOAT, nnod);
    mesh->vvbar = ARRAY(CHFLOAT, nnod);
    mesh->xv05 = ARRAY(CHFLOAT, nnod);
    mesh->yv05 = ARRAY(CHFLOAT, nnod);
    mesh->pre05 = ARRAY(CHFLOAT, nel);
    mesh->rho05 = ARRAY(CHFLOAT, nel);
    mesh->en05 = ARRAY(CHFLOAT, nel);
    mesh->volel05 = ARRAY(CHFLOAT, nel);
    mesh->divint = ARRAY(CHFLOAT, nel);
    mesh->divvel = ARRAY(CHFLOAT, nel);

    SETARRAY2D(mesh->nint, CHFLOAT, nel, 4);
    SETARRAY2D(mesh->elwtc, CHFLOAT, nel, 4);
    SETARRAY2D(mesh->dndx, CHFLOAT, nel, 4);
    SETARRAY2D(mesh->dndy, CHFLOAT, nel, 4);
    SETARRAY2D(mesh->pdndx, CHFLOAT, nel, 4);
    SETARRAY2D(mesh->pdndy, CHFLOAT, nel, 4);

    return 0;
}

int main() {
    char meshfile[CHSTRLEN] = "/home/nick/Programming/Chocobo_Runs/ChocoboF/sod/mesh.chc";
    char eosfile[CHSTRLEN] = "/home/nick/Programming/Chocobo_Runs/ChocoboF/sod/eos.dat";
    CHINT nel, nnod, nreg, nmat;

    struct MeshT mesh;
    struct InputT input = {0.75, 0.5, 0.0001, 0.0001, 1.02, 0.0, 0.205, 0.2, 0, 0};

    // Get mesh size
    get_mesh_size(meshfile, &nel, &nnod, &nreg, &nmat);

    // Allocate storage
    set_store(&mesh, nel, nnod, nreg, nmat);

    // Read mesh
    read_mesh(meshfile, nel, nnod, mesh.nodelist, mesh.znodbound, mesh.region,
              mesh.material, mesh.regiontocell, mesh.xv, mesh.yv);

    // Set initial conditions
    set_initial_conditions(eosfile, nel, nmat, mesh.material, mesh.pre,
                           mesh.rho, mesh.gamma);

    // Initial volume
    calculate_volume(nel, nnod, mesh.xv, mesh.yv, mesh.nodelist, mesh.volel,
                     mesh.area);

    // Initial mass
    calculate_mass(nel, mesh.volel, mesh.rho, mesh.massel);

    // Initial energy
    CELL_LOOP {
        int mat = mesh.material[cell];
        mesh.en[cell] = mesh.pre[cell]/((mesh.gamma[mat] - 1.0)*mesh.rho[cell]);
    }

    // Calculate total energy
    CHFLOAT total_energy, total_ie, total_ke;
    get_total_energy(nel, mesh.en, mesh.rho, mesh.uv, mesh.vv, mesh.massel,
                          mesh.elwtc, mesh.nodelist,
                          &total_energy, &total_ke, &total_ie);

    printf("\nInitial energy::\nTotal: %f\nInternal: %f\nKinetic: %f\n\n",
           total_energy, total_ie, total_ke);

    // Main loop
    CHFLOAT t = input.t0;
    CHFLOAT dt = input.dtinit;
    CHINT step = 0;
    CHINT dtcontrol;
    while (t < input.tend) {
        step += 1;
        printf("Step: %d/%d", step, input.stepcount);

        // STUFF
        // Finite Elements
        calculate_finite_elements(nel, mesh.xv, mesh.yv, mesh.nodelist, mesh.nint,
                                  mesh.dndx, mesh.dndy, mesh.pdndx, mesh.pdndy,
                                  mesh.elwtc);

        // Divergence of V
        calculate_div_v(nel, mesh.uv, mesh.vv, mesh.pdndx, mesh.pdndy,
                        mesh.nodelist, mesh.divvel);

        // Soundspeed
        calculate_soundspeed(nel, mesh.pre, mesh.rho, mesh.material, mesh.gamma,
                            mesh.cc);

        // Artificial viscosity
        calculate_q(nel, mesh.rho, mesh.cc, mesh.divvel, mesh.area, input.cq,
                    input.cl, mesh.qq);

        // Timestep
        get_dt(nel, mesh.rho, mesh.area, mesh.cc, mesh.qq, t, input.t0,
               input.dtinit, input.dtmax, input.growth, &dt, &dtcontrol);

        t += dt;
        printf(", time: %f, dt: %f, Element: %d\n", t, dt, dtcontrol);

        // Half timestep
        CHFLOAT dt05 = dt/2.0;

        // Half timestep nodal positions
        move_nodes(nnod, dt05, mesh.xv, mesh.yv, mesh.uv, mesh.vv, mesh.xv05,
                   mesh.yv05);

        // Finite Elements - Half timestep
        calculate_finite_elements(nel, mesh.xv05, mesh.yv05, mesh.nodelist,
                                  mesh.nint, mesh.dndx, mesh.dndy, mesh.pdndx,
                                  mesh.pdndy, mesh.elwtc);

        // Store volume
        CELL_LOOP {
            mesh.volelold[cell] = mesh.volel[cell];
        }

        // Half timestep volume
        calculate_volume(nel, nnod, mesh.xv05, mesh.yv05, mesh.nodelist,
                         mesh.volel05, mesh.area);

        // Half timestep density
        calculate_density(nel, mesh.massel, mesh.volel05, mesh.rho05);

        // Integral of divergence of V
        calculate_int_divv(input.zintdivvol, nel, nnod, dt05, mesh.volel,
                           mesh.volelold, mesh.uv, mesh.vv, mesh.dndx,
                           mesh.dndy, mesh.nodelist, mesh.divint);

        // Half timestep energy
        calculate_energy(nel, dt05, mesh.pre, mesh.qq, mesh.massel, mesh.en,
                         mesh.divint, mesh.en05);

        // Half timestep pressure from EoS
        perfect_gas(nel, mesh.en05, mesh.rho05, mesh.material, mesh.gamma,
                    mesh.pre05);

        // Store old velocities
        NODE_LOOP {
            mesh.uvold[node] = mesh.uv[node];
            mesh.vvold[node] = mesh.vv[node];
        }

        // Momentum calculation
        momentum_calculation(nel, nnod, dt, 0, 0, 0.0,
                             mesh.uv, mesh.vv, mesh.xv,
                             mesh.yv, mesh.rho,
                             mesh.pre, mesh.area,
                             mesh.cc, mesh.qq,
                             mesh.nint, mesh.dndx,
                             mesh.dndy, mesh.nodelist,
                             mesh.znodbound, mesh.uv,
                             mesh.vv);

        // Average velocities
        NODE_LOOP {
            mesh.uvbar[node] = (mesh.uvold[node] + mesh.uv[node])/2.0;
            mesh.vvbar[node] = (mesh.vvold[node] + mesh.vv[node])/2.0;
        }

        // Full timestep node positions
        move_nodes(nnod, dt, mesh.xv, mesh.yv, mesh.uv, mesh.vv, mesh.xv,
                   mesh.yv);

        // Finite Elements
        calculate_finite_elements(nel, mesh.xv, mesh.yv, mesh.nodelist, mesh.nint,
                                  mesh.dndx, mesh.dndy, mesh.pdndx, mesh.pdndy,
                                  mesh.elwtc);

        // Full timestep volume
        calculate_volume(nel, nnod, mesh.xv, mesh.yv, mesh.nodelist,
                         mesh.volel, mesh.area);

        // Full timestep density
        calculate_density(nel, mesh.massel, mesh.volel, mesh.rho);

        // Integral of divergence of V
        calculate_int_divv(input.zintdivvol, nel, nnod, dt, mesh.volel,
                           mesh.volelold, mesh.uvbar, mesh.vvbar, mesh.dndx,
                           mesh.dndy, mesh.nodelist, mesh.divint);

        // Full timestep energy
        calculate_energy(nel, dt, mesh.pre05, mesh.qq, mesh.massel, mesh.en,
                         mesh.divint, mesh.en);

        // Full timestep pressure from EoS
        perfect_gas(nel, mesh.en, mesh.rho, mesh.material, mesh.gamma, mesh.pre);

        if (step >= input.stepcount && input.stepcount > 0) break;
    }
}
