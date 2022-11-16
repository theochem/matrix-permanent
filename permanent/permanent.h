#ifndef PERMANENT_H
#define PERMANENT_H


#include <stdint.h>


double opt(const int64_t, const int64_t, const double *);


double combinatoric(const int64_t, const int64_t, const double *);


double glynn(const int64_t, const int64_t, const double *);


double ryser(const int64_t, const int64_t, const double *);


double glynn_rectangle(const int64_t, const int64_t, const double *);


double ryser_rectangle(const int64_t, const int64_t, const double *);


#endif /* PERMANENT_H */
