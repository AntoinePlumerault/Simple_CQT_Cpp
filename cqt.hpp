#ifndef CQT_HPP
#define CQT_HPP
#include "fft.hpp"
#include <vector>

class CQT {
public:
	CQT(size_t n, size_t n_bins, double f_min, double f_max, double f_e);
	void cqt(Complex & x);

private:
	size_t m_N;
	size_t m_K;
	size_t* m_Delta_omega_cqt;
	std::vector<double*> m_kernels;

};
#endif
