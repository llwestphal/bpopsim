#include "pop_sim_header.h"
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
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost;
using namespace boost::program_options;

mt19937 rng(time(NULL));

double returnExp(double a)
{
	exponential_distribution<> exp_dist(a);
	variate_generator<mt19937&, exponential_distribution <> > next_value(rng, exp_dist);
	return next_value(); //Not sure why I need to do this...
}

uint64_t returnBin(const uint64_t n, const double p)
{
	binomial_distribution<> my_binomial(n,p);
	variate_generator<mt19937&, binomial_distribution<> > next_value(rng, my_binomial);
	return next_value();
}

// NOTE: The BOOST Poisson random number generator DOES NOT WORK for large means!!!
double returnPoisson(const double n, const double p)
{
	poisson_distribution<> my_poisson(n * p);
	variate_generator<mt19937&, poisson_distribution<> > next_value(rng, my_poisson);
	cout << "Poisson with mean " << n << " * " << p << " = " << next_value() << endl;
	return next_value();
}
int cline(int argc, char* argv[]) {
// setup and parse configuration options:
	options_description cmdline_options("Allowed options");
        cmdline_options.add_options()
        ("help,h", "produce this help message")
        ("bam,b", value<string>(), "bam file containing sequences to be aligned")
        ("fasta,f", value<string>(), "FASTA file of reference sequence")
        ("output,o", value<string>(), "output directory")
        ("readfile,r", value<vector<string> >(), "names of readfiles (no extension)");

        variables_map options;
        store(parse_command_line(argc, argv, cmdline_options), options);
        notify(options);

        // make sure that the config options are good:
        if(options.count("help")
                 || !options.count("bam")
                 || !options.count("fasta")
                 || !options.count("output")
                 || !options.count("readfile")) {
                cout << "Usage: error_count --bam <sequences.bam> --fasta <reference.fasta> --output <path> --readfile <filename>" << endl;
                cout << cmdline_options << endl;
                return -1;
        }
        return 0;
}
