#ifndef common_h
#define common_h

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <iomanip> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <stdint.h>
#include <time.h>
#include "tree.hh"
#include "tree_util.hh"
#include "icsilog.h"
#include <functional>
#include <list>
#include <utility>
#include <map>

#include <boost/program_options.hpp>

// Global variable for keeping track of verbosity
extern bool g_verbose;
extern bool g_ro_only;


#endif
