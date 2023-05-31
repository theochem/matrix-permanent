#ifndef PERMANENT_H
#define PERMANENT_H


#include <stdint.h>


double opt(const int64_t, const int64_t, double *const);


double combinatoric(const int64_t, const int64_t, double *const);


double glynn(const int64_t, const int64_t, double *const);


double ryser(const int64_t, const int64_t, double *const);


double glynn_rectangle(const int64_t, const int64_t, double *const);


double ryser_rectangle(const int64_t, const int64_t, double *const);


#endif /* PERMANENT_H */
