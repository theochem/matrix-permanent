#include <complex>

#include "permanent.h"


template<> double combinatoric<double>(const std::size_t, const std::size_t, double *const);
template<> double combinatoric_rectangular<double>(const std::size_t, const std::size_t, double *const);
template<> double glynn<double>(const std::size_t, const std::size_t, double *const);
template<> double glynn_rectangular<double>(const std::size_t, const std::size_t, double *const);
template<> double ryser<double>(const std::size_t, const std::size_t, double *const);
template<> double ryser_rectangular<double>(const std::size_t, const std::size_t, double *const);
template<> double opt<double>(const std::size_t, const std::size_t, double *const);

template<> std::complex<double> combinatoric<std::complex<double>>(const std::size_t, const std::size_t, std::complex<double> *const);
template<> std::complex<double> combinatoric_rectangular<std::complex<double>>(const std::size_t, const std::size_t, std::complex<double> *const);
template<> std::complex<double> glynn<std::complex<double>>(const std::size_t, const std::size_t, std::complex<double> *const);
template<> std::complex<double> glynn_rectangular<std::complex<double>>(const std::size_t, const std::size_t, std::complex<double> *const);
template<> std::complex<double> ryser<std::complex<double>>(const std::size_t, const std::size_t, std::complex<double> *const);
template<> std::complex<double> ryser_rectangular<std::complex<double>>(const std::size_t, const std::size_t, std::complex<double> *const);
template<> std::complex<double> opt<std::complex<double>>(const std::size_t, const std::size_t, std::complex<double> *const);
