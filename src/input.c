#include "ini.h"
#include "chocobo_core.h"
#include "data.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

static int handler(void* input, const char* section, const char* name,
                   const char* value)
{
    InputTT* config = (InputTT*)input;

    if (MATCH("Q","CL")) {
        config->cl = atof(value);
    }
    else if (MATCH("Q","CQ")) {
        config->cq = atof(value);
    }

    else if (MATCH("Control","t0")) {
        config->t0 = atof(value);
    }
    else if (MATCH("Control","tend")) {
        config->tend = atof(value);
    }
    else if (MATCH("Control","dtinit")) {
        config->dtinit = atof(value);
    }
    else if (MATCH("Control","dtmax")) {
        config->dtmax = atof(value);
    }
    else if (MATCH("Control","growth")) {
        config->growth = atof(value);
    }

    else if (MATCH("Debug","debug_step_count")) {
        config->stepcount = atof(value);
    }

    else if (MATCH("Material","material_numbers")) {
        // Currently ignored
    }

    else if (strcasecmp(section, "SphSod Mesh") == 0) {
        // Ignore the hard-coded mesh input, as we don't yet support that.
    }

    else {
        printf("Unexpected variable %s in section %s\n", name, section);
        exit(-1);
    }

    return 0;

}

CHRETURN read_input(char inputfile[CHSTRLEN], InputTT* input)
{
    ini_parse(inputfile, handler, input);
}
