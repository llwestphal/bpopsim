#ifndef POP_SIM_HEADER_H
#define POP_SIM_HEADER_H

/*pop_sim_header.h*/
#include <stdint.h>

double returnExp(double a);
uint64_t returnBin(const uint64_t n, const double p);

// NOTE: The BOOST Poisson random number generator DOES NOT WORK for large means!!!

struct Individual {
	double w;
	double n;
	char color;
};

double returnPoisson(const double n, const double p);

int cline(int argc, char* argv[]);



#endif
