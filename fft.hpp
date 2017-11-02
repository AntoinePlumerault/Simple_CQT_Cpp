#ifndef FFT_HPP
#define FFT_HPP

#include <complex> 
#include <valarray>

#define PI 3.14159265359

typedef std::complex<double> Complex;

void separate(Complex* X, const size_t N);
void fft(Complex* X, const size_t N);

#endif 

