#ifndef RAM_H

#define RAM_H

using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <deque>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>
#include <assert.h>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash/hash.hpp>
#include <boost/algorithm/string.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
     
#include <getopt.h>
     
#define UMILLION 1000000ULL
#ifdef WINDOWS
	void gettimeofday( struct timeval* p, void* );
#endif

#include "genotype.h"
#include "io.h"
#include "stringfn.h"
#include "mathfn.h"
#include "data.h"
#include "fileio.h"
#include "functions.h"

#endif
