#ifndef CQT_HPP
#define CQT_HPP

#include <vector>
#include <iostream>
#include <complex>

#include "fft.hpp"

#define PI 3.14159265359

typedef std::complex<double> Complex;

class CQT
{
public:
	CQT();
	CQT(size_t n, size_t n_bins, double f_min, double f_max, double f_e);
	Complex * transform(Complex * x) const;

private:
	size_t m_N;
	size_t m_K;
	size_t* m_Delta_cqt;
	double* m_omega_cqt;
	double m_omega_e;
	std::vector<double*> m_kernels;
};

#endif
