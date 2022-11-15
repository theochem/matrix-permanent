#ifndef PERMANENT_H
#define PERMANENT_H


#include <stdint.h>


#ifdef TUNING_FILE
/* Include the tuning file. */
#include "tuning.h"
#else
/* Set default tuning parameters. */
#define PARAM_1 3.1415926
#define PARAM_2 8192
#endif


double opt(const int64_t, const int64_t, const double *);


double combinatoric(const int64_t, const int64_t, const double *);


double glynn(const int64_t, const int64_t, const double *);


double ryser(const int64_t, const int64_t, const double *);


#endif /* PERMANENT_H */
