#include<iostream>
#include<iomanip>
#include<cmath>

// x^(y) % p
double modexp(double x, double y, double p)
{
	double res = 1;

	//Update x if it is greater than or equal to p
	x = std::fmod(x, p);

	while (y > 0)
	{
		// If y is odd, multiply x with result
		if (((int)y & 1) == 1)
			res = std::fmod(res * x, p);
		// y must be even now
		y = (int)y >> 1;
		x = std::fmod(x * x, p);
	}
	return res;
}

double expm(double p, double ak)
{
	int i, j;
	double p1, pt, r;
#define ntp 25
	static double tp[ntp];
	static int tp1 = 0;

	if (tp1 == 0) {
		tp1 = 1;
		tp[0] = 1.;

		for (i = 1; i < ntp; i++) tp[i] = 2. * tp[i - 1];
	}

	if (ak == 1.) return 0.;

	for (i = 0; i < ntp; i++) if (tp[i] > p) break;

	pt = tp[i - 1];
	p1 = p;
	r = 1.;

	for (j = 1; j <= i; j++) {
		if (p1 >= pt) {
			r = 16. * r;
			r = r - (int)(r / ak) * ak;
			p1 = p1 - pt;
		}
		pt = 0.5 * pt;
		if (pt >= 1.) {
			r = r * r;
			r = r - (int)(r / ak) * ak;
		}
	}

	return r;
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
	double pid =  4 * getS_j(n, 1) - 2 * getS_j(n, 4) - getS_j(n, 5) - getS_j(n, 6);
	double fractionalpart = std::modf(pid, &fractionalpart);
	pid = pid - (int)pid + 1.0;
	return pid;
}

int main()
{
	//calculate the 0th digit of pi
	std::cout << std::setprecision(15) << getnthpi(0) << std::endl;
	//calculate the 1000000th digit of pi
	std::cout << std::setprecision(15) << getnthpi(1000000) << std::endl;
	return 0;
}

