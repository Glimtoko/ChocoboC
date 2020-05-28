#include "chocobo_core.h"

CHRETURN output(char baseloc[CHSTRLEN], CHINT step, CHINT nel, CHINT nnod,
                CHFLOAT time, CHFLOAT xv[nnod], CHFLOAT yv[nnod],
                CHFLOAT uv[nnod], CHFLOAT vv[nnod], CHINT **nodelist,
                CHFLOAT rho[nel], CHFLOAT pre[nel], CHFLOAT en[nel],
                CHFLOAT volume[nel], CHINT *filenum);
