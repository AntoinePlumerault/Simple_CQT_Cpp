#include "fft.hpp"
#include <math.h>

void separate(Complex* X, const size_t N) {
	Complex* temp = new Complex[N / 2];
	for (size_t i = 0; i < N / 2; ++i)
	{
		temp[i] = X[2 * i + 1];
		X[i] = X[2 * i];
	}
	for (size_t i = 0; i < N / 2; i++)
	{
		X[i + N / 2] = temp[i];
	}
	delete temp;
}

void fft(Complex* X, const size_t N) {
	if (N < 2)
	{
		// Do nothing
	}
	else
	{
		separate(X, N);     //even numbers to the left, odd numbers to the left 
		fft(X, N / 2); 
		fft(X + N / 2, N / 2);

		for (size_t i = 0; i < N / 2; ++i)
		{
			Complex even = X[k];   
			Complex odd  = X[k + N / 2];
			Complex w    = std::exp(Complex(0, -2.*PI*i / N));
			X[k]         = even + w * odd;
			X[k + N / 2] = even - w * odd;
		}
	}
}
