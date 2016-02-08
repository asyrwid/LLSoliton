#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H
#include <complex>
#include <vector>
#include <iostream>
#include <omp.h>
using std::vector;
using std::cout;

class decomposition
{


public:

    decomposition();

    vector<vector<double>> MatrixOfZeros(int rows, int columns);

    void Coefficients_CramerMethod(int NumOfParticlesBeforeMeasurement,
                                   vector<double> MeasuredPositions,
                                   int NumOfStates,
                                   double L);
};

#endif // DECOMPOSITION_H
