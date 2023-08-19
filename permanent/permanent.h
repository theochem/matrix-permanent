#ifndef PERMANENT_H
#define PERMANENT_H

#include <stdlib.h>


double opt(const size_t, const size_t, double *const);


double combinatoric(const size_t, const size_t, double *const);


double glynn(const size_t, const size_t, double *const);


double ryser(const size_t, const size_t, double *const);


double combinatoric_rectangle(const size_t, const size_t, double *const);


double glynn_rectangle(const size_t, const size_t, double *const);


double ryser_rectangle(const size_t, const size_t, double *const);


#endif /* PERMANENT_H */
