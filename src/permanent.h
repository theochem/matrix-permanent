#ifndef PERMANENT_PERMANENT_H
#define PERMANENT_PERMANENT_H


#include <stdlib.h>


double opt(const size_t m, const size_t n, double *const ptr);


double combinatoric(const size_t m, const size_t n, double *const ptr);


double combinatoric_rectangular(const size_t m, const size_t n, double *const ptr);


double glynn(const size_t m, const size_t n, double *const ptr);


double glynn_rectangular(const size_t m, const size_t n, double *const ptr);


double ryser(const size_t m, const size_t n, double *const ptr);


double ryser_rectangular(const size_t m, const size_t n, double *const ptr);


#endif /* PERMANENT_PERMANENT_H */
