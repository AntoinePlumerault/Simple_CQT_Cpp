#include "cqt.hpp"

CQT::CQT(size_t N, size_t K, double f_min, double f_max, double f_e) {
	m_N = N;
	size_t n_octave = round(log(f_max / f_min) / log(2));
	size_t n_bins = K / n_octave;
	double Q = 1.0 / (pow(2.0, 1.0 / n_bins));
	double sigma_cqt;
	double omega_cqt;
	m_Delta_omega_cqt = new size_t[K];
	for (size_t k_cq = 0; k_cq < K; ++k_cq)
	{
		omega_cqt = 2.0*PI*f_min * (pow(2.0, k_cq / n_bins));
		sigma_cqt = N*omega_cqt / (PI * Q * f_e);
		m_Delta_omega_cqt[k_cq] = ceil(6.0 * sigma_cqt);
		m_kernels[k_cq] = new double[m_Delta_omega_cqt[k_cq]];
		
		for (size_t k = 0; k < m_Delta_omega_cqt[k_cq]; ++k)
		{
			m_kernels[k_cq][k] = exp(-1.0 / 2.0 * pow((k - m_Delta_omega_cqt[k_cq] / 2) / (sigma_cqt), 2.0));
		}
	}

}
void CQT::cqt(Complex & x) {
	fft(x, m_N);
	for (size_t j = 0; j < N; ++j)
	{

	}
}