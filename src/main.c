#include "chocobo_core.h"
#include "hydro.h"
#include "data.h"
#include "mesh.h"

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

//     mesh->nodelist = ARRAY(CHINT, nel*4);
    SETARRAY2D(mesh->nodelist, CHINT, nel, 4);
    mesh->region = ARRAY(CHINT, nel);
    mesh->material = ARRAY(CHINT, nel);
    mesh->regiontocell = ARRAY(CHINT, nel*2);

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

    mesh->nint = ARRAY(CHFLOAT, nel*4);
    mesh->elwtc = ARRAY(CHFLOAT, nel*4);
    mesh->dndx = ARRAY(CHFLOAT, nel*4);
    mesh->dndy = ARRAY(CHFLOAT, nel*4);
    mesh->pdndx = ARRAY(CHFLOAT, nel*4);
    mesh->pdndy = ARRAY(CHFLOAT, nel*4);

    return 0;
}

int main() {
    char meshfile[CHSTRLEN] = "/home/nick/Programming/Chocobo_Runs/ChocoboF/sod_large/mesh.chc";
    CHINT nel, nnod, nreg, nmat;

    struct MeshT mesh;
    struct InputT input = {0.75, 0.5, 0.0001, 0.0001, 1.02, 0.0, 0.205, 0.2, 0, 1};

    // Get mesh size
    get_mesh_size(meshfile, &nel, &nnod, &nreg, &nmat);

    // Allocate storage
    set_store(&mesh, nel, nnod, nreg, nmat);

    read_mesh(meshfile, nel, nnod, mesh.nodelist, mesh.region, mesh.material,
              mesh.regiontocell, mesh.xv, mesh.yv);

    // Main loop
    CHFLOAT t = input.t0;
    CHINT step = 0;
    while (t < input.tend) {
        step += 1;
        printf("Step %d/%d\n", step, input.stepcount);

        // STUFF

        if (step >= input.stepcount && input.stepcount > 0) break;
    }
}
