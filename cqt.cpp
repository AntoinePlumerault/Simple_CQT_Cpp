#include "cqt.hpp"



CQT::CQT()
{

}

CQT::CQT(size_t N, size_t K, double f_min, double f_max, double f_e) : CQT()
{
	size_t n_octave = (size_t)round(log(f_max / f_min) / log(2));
	size_t n_bins = K / n_octave;
	double Q = 1.0 / (pow(2.0, 1.0 / n_bins) - 1.0);
	double sigma_cqt;
	m_N = N;
	m_K = K;
	m_kernels.resize(m_K);
	m_omega_e = 2.0*PI*f_e;
	m_omega_cqt = new double[m_K];
	m_Delta_cqt = new size_t[m_K];
	for (size_t k_cq = 0; k_cq < K; ++k_cq)
	{
		m_omega_cqt[k_cq] = 2.0*PI*f_min * (pow(2.0, (double)k_cq / (double)n_bins));
		sigma_cqt = (double)N*m_omega_cqt[k_cq] / (PI * Q * 2.0*PI*f_e);
		m_Delta_cqt[k_cq] = (size_t)ceil(6.0 * sigma_cqt / 2.0 );
		
		//std::cout << N << "   " << Q << "   " << n_bins << "   " << k_cq << "   " << sigma_cqt << "   " << m_omega_cqt[k_cq]  << "   " << m_Delta_cqt[k_cq] << std::endl;
		//std::cout << 2.0*m_N*m_omega_cqt[k_cq] / m_omega_e << std::endl;
		
		double * kernel;
		kernel = new double[2 * m_Delta_cqt[k_cq] + 1];
		m_kernels[k_cq] = kernel;
		
		//std::cout << 222222222222 << std::endl;

		for (size_t k = 0; k < 2*m_Delta_cqt[k_cq] + 1; ++k)
		{
			m_kernels[k_cq][k] = exp(- 1.0 / 2.0 * pow((k - m_Delta_cqt[k_cq]) / (sigma_cqt), 2.0));
		}
	}
}


Complex* CQT::transform(Complex * x) const
{
	fft(x, 2*m_N);

	Complex * Qx;
	Qx = new Complex[m_K];

	for (size_t k_cq = 0; k_cq < m_K; ++k_cq)
	{
		Qx[k_cq] = 0.0;
		//std::cout << -(int)m_Delta_cqt[k_cq] << "   " << (int)m_Delta_cqt[k_cq] << std::endl;
		for (int k = -(int)m_Delta_cqt[k_cq]; k < (int)m_Delta_cqt[k_cq] + 1; ++k)
		{
			//std::cout << "bbbbbbbbb" << std::endl;
			Qx[k_cq] += m_kernels[k_cq][k + m_Delta_cqt[k_cq]] * x[(size_t)(round(k + 2.0*m_N*m_omega_cqt[k_cq] / m_omega_e))];
		}
		
		//Qx[k_cq] = Qx[k_cq] / (double)m_N;
	}
	return Qx;
}