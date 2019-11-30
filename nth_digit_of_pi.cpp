/*
https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula
https://www.davidhbailey.com/dhbpapers/bbp-alg.pdf
*/

#include<iostream>
#include<iomanip>
#include<cmath>

// calculate x^(y) % p, https://www.geeksforgeeks.org/modular-exponentiation-power-in-modular-arithmetic/
double modexp(double x, double y, double p)
{
	double remainder = 1;

	//Update x if it is greater than or equal to p
	x = std::fmod(x, p);

	while (y > 0)
	{
		// If y is odd, multiply x with result
		if (((int)y & 1) == 1)
			remainder = std::fmod(remainder * x, p);
		// y must be even now
		y = (int)y >> 1;
		x = std::fmod(x * x, p);
	}
	return remainder;
}

double getS_j(int n, int j)
{
	double sum = 0;

	for (int k = 0; k < n; k++)
	{
		sum = sum + modexp(16, n - k, (8 * k + j)) / (8 * k + j);
		//sum = sum + expm(n - k, (8 * k + j)) / (8 * k + j);
		sum = sum - (int)sum;
	}
	
	const double eps = 1e-17;

	for (long k = n; k < n + 100; k++)
	{
		double quotient = std::pow(16.0, (double)n - k) / (8 * k + j);
		if (quotient < eps)
			break;
		sum = sum + quotient;
		sum = sum - (int)sum;
	}
	
	return sum;
}

double getnthpi(long n)
{
	double fractionalpart = std::modf(4 * getS_j(n, 1) - 2 * getS_j(n, 4) - getS_j(n, 5) - getS_j(n, 6), &fractionalpart);
	fractionalpart += 1.0;
	return fractionalpart;
}

int main()
{
	/* last digit displayed is not accurate */
	std::cout << std::setprecision(13) << getnthpi(0) << std::endl; //calculate the 0th digit of pi
	std::cout << std::setprecision(13) << getnthpi(1000000) << std::endl; //calculate the 1000000th digit of pi
	return 0;
}

