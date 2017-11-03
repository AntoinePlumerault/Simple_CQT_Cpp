#include <math.h>
#include "cqt.hpp"

CQT::CQT(size_t N, size_t K, double f_min, double f_max, double f_e) {
	m_N = N;
	size_t n_octave = round(log(f_max / f_min) / log(2));
	size_t n_bins = K / n_octave;
	double Q = 1.0 / (pow(2.0, 1.0 / n_bins));
	double sigma_cqt;
	m_omega_e = 2.0*PI*f_e;
	m_omega_cqt = new double[K];
	m_Delta_cqt = new size_t[K];
	for (size_t k_cq = 0; k_cq < K; ++k_cq)
	{
		m_omega_cqt[k_cq] = 2.0*PI*f_min * (pow(2.0, k_cq / n_bins));
		sigma_cqt = N*m_omega_cqt[k_cq] / (PI * Q * f_e);
		m_Delta_cqt[k_cq] = ceil(6.0 * sigma_cqt / 2.0);
		m_kernels[k_cq] = new double[m_Delta_cqt[k_cq]];
		
		for (size_t k = 0; k < 2*m_Delta_cqt[k_cq] + 1; ++k)
		{
			m_kernels[k_cq][k] = exp(-1.0 / 2.0 * pow((k - m_Delta_cqt[k_cq]) / (sigma_cqt), 2.0));
		}
	}

}
double* CQT::cqt(Complex* x) {
	fft(x, m_N);
	double* X = new double[m_K];
	for (size_t k_cq = 0; k_cq < m_K; ++k_cq)
	{
		for (int k = -(int)m_Delta_cqt[k_cq]; k < m_Delta_cqt[k_cq] + 1; ++k)
		{
			X[k_cq] += (m_kernels[k_cq][k + m_Delta_cqt[k_cq]] * x[(size_t)(round(k + 2.0*m_N*m_omega_cqt[k_cq] / m_omega_e))]).real();
		}
		X[k_cq] /= m_N;
	}
	return X;
}