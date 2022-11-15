#include <stdint.h>
#include <stdio.h>

#include "permanent.h"


int main()
{
    /* Compute the tuning parameters */
    double param_1 = 3.14159;
    int64_t param_2 = 8192;
    /* ...
     * ...
     * ... */

    /* Print a header file to stdout with constants defined as macros */
    printf("#ifndef PERMANENT_TUNING_H\n");
    printf("#define PERMANENT_TUNING_H\n");
    printf("\n\n");
    printf("#define PARAM_1 %.9le\n", param_1);
    printf("#define PARAM_2 %ld\n", param_2);
    printf("\n\n");
    printf("#endif /* PERMANENT_TUNING_H */\n");

    /* Exit successfully */
    return 0;
}
