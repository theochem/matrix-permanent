/* Copyright 2034 QC-Devs (GPLv3) */

#include "permanent.h"

#include <complex>

template double combinatoric<double>(const std::size_t, const std::size_t,
                                     const double *);

template double combinatoric_rectangular<double>(const std::size_t,
                                                 const std::size_t,
                                                 const double *);

template double glynn<double>(const std::size_t, const std::size_t,
                              const double *);

template double glynn_rectangular<double>(const std::size_t, const std::size_t,
                                          const double *);

template double ryser<double>(const std::size_t, const std::size_t,
                              const double *);

template double ryser_rectangular<double>(const std::size_t, const std::size_t,
                                          const double *);

template double opt<double>(const std::size_t, const std::size_t,
                            const double *);

template std::complex<double> combinatoric<std::complex<double>>(
    const std::size_t, const std::size_t, const std::complex<double> *);

template std::complex<double> combinatoric_rectangular<std::complex<double>>(
    const std::size_t, const std::size_t, const std::complex<double> *);

template std::complex<double> glynn<std::complex<double>>(
    const std::size_t, const std::size_t, const std::complex<double> *);

template std::complex<double> glynn_rectangular<std::complex<double>>(
    const std::size_t, const std::size_t, const std::complex<double> *);

template std::complex<double> ryser<std::complex<double>>(
    const std::size_t, const std::size_t, const std::complex<double> *);

template std::complex<double> ryser_rectangular<std::complex<double>>(
    const std::size_t, const std::size_t, const std::complex<double> *);

template std::complex<double> opt<std::complex<double>>(
    const std::size_t, const std::size_t, const std::complex<double> *);
