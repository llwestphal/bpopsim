#include "random.h"

using namespace std;
using namespace boost;

	double RandomGenerator::expDist(double f)

	{
	}
	uint64_t RandomGenerator::binDist(const uint64_t d, const double e)
	{
	}
	double RandomGenerator::poisDist(const double a, const double b)
	{
	}


//mt19937 rng(time(NULL));
//
//double returnExp(double a)
//{
//	exponential_distribution<> exp_dist(a);
//	variate_generator<mt19937&, exponential_distribution <> > next_value(rng, exp_dist);
//	return next_value(); //Not sure why I need to do this...
//}
//
//uint64_t returnBin(const uint64_t n, const double p)
//{
//	binomial_distribution<> my_binomial(n,p);
//	variate_generator<mt19937&, binomial_distribution<> > next_value(rng, my_binomial);
//	return next_value();
//}
//
// NOTE: The BOOST Poisson random number generator DOES NOT WORK for large means!!!
//double returnPoisson(const double n, const double p)
//{
//	poisson_distribution<> my_poisson(n * p);
//	variate_generator<mt19937&, poisson_distribution<> > next_value(rng, my_poisson);
//	cout << "Poisson with mean " << n << " * " << p << " = " << next_value() << endl;
//	return next_value();
//}
