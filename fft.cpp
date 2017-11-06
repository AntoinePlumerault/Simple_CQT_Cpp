#include "fft.hpp"

void separate(Complex * X, size_t N) 
{
	Complex * b;
	b = new Complex[N/2];

	for (size_t i = 0; i < N / 2; ++i)
	{
		b[i] = X[2 * i + 1];
		X[i] = X[2 * i];
	}
	for (size_t i = 0; i < N / 2; ++i)
	{
		X[i + N / 2] = b[i];
	}
	delete [] b;
}

void fft(Complex * X, size_t N) 
{
	if (N < 2)
	{
		// Do nothing
	}
	else
	{
		separate(X, N);     //even numbers to the left, odd numbers to the left 
		fft(X,         N / 2); 
		fft(X + N / 2, N / 2);

		for (size_t k = 0; k < N / 2; ++k)
		{
			Complex even = X[k];   
			Complex odd  = X[k + N / 2];
			Complex w    = std::exp(Complex(0.0, -2.*PI*k / N));
			X[k]         = even + w * odd;
			X[k + N / 2] = even - w * odd;
		}
	}
}
