#ifndef UTILS
#define UTILS
#include <iterator>
#include <math.h>
#include <random>
#include <omp.h>
#include <iostream>
#include <fstream>

# define M_PI  3.14159265358979323846  /* pi definition */

namespace utils
{


double GetRandom(double min, double max)
{// Returns a random double between min and max
    return ((double) rand()*(max-min)/(double)RAND_MAX + min);
}



}
#endif // UTILS

