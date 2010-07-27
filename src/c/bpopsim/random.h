#ifndef RANDOM_
#define RANDOM_

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#if !defined(__SUNPRO_CC) || (__SUNPRO_CC > 0x530)
#include <boost/generator_iterator.hpp>
#endif
#include <boost/math/special_functions/factorials.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/random/poisson_distribution.hpp>

class RandomGenerator

{

	public:
	RandomGenerator();
	double expDist(double);
	uint64_t binDist(const uint64_t, const double);
	double poisDist(const double, const double);
};

#endif
