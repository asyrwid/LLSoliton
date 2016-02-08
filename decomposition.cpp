#include "decomposition.h"
#include "permutator.h"
#include "utils.h"
#include <iterator>
#include <math.h>
#include <random>
#include <omp.h>
#include <iostream>
#include <fstream>

typedef std::complex<double> Complex;

decomposition::decomposition()
{
}

vector<vector<double>> decomposition::MatrixOfZeros(int rows, int columns)
{
    vector<double> tmp;
    for(int j = 0; j< columns; j++)
    {
        tmp.push_back(0.);
    }

    vector<vector<double>> A;
    for(int i = 0; i < rows; i ++)
    {
        A.push_back(tmp);
    }

    return A;
}


void decomposition::Coefficients_CramerMethod(int NumOfParticlesBeforeMeasurement,
                                              vector<double> MeasuredPositions,
                                              int NumOfStates,
                                              double L)
{
    // Collections of positions preparation
    //==============================================================================================
    vector<vector<double>> MatrixX0, MatrixX1;
    vector<double> X1, X0;
    int NumOfParticlesAfterMeasurement = NumOfParticlesBeforeMeasurement - MeasuredPositions.size();
    X0 = MeasuredPositions;
    double epsilon = .1;

    for(int u = 0; u < NumOfParticlesAfterMeasurement - 1; u++)
    {
        double delta = utils::GetRandom( -epsilon,  epsilon);
        X0.push_back(MeasuredPositions[0] + delta);
        X1.push_back(MeasuredPositions[0] + delta);
    }

    for(int y = 0; y < NumOfStates; y++)
    {
        vector<double> tmp0, tmp1;
        tmp0 = X0;
        tmp1 = X1;
        tmp0.push_back(  1.* y/(NumOfStates - 1)  );
        tmp1.push_back(  1.* y/(NumOfStates - 1)  );
        MatrixX0.push_back( tmp0 );
        MatrixX1.push_back( tmp1 );
    }
    //==============================================================================================

    //We can choose arbitrary phase factor. In the case of initial state I choose Exp(I * .5)
    Complex Phase_Factor;
    Phase_Factor = exp( Complex(0., .5 ) );
    //==============================================================================================

    //Matrices preparation - we have 2 x NumOfStates parameters - real and imaginary part
    vector<vector<vector<double>>> ReA, ImA;

    for(int i = 0; i < NumOfStates; i++)
    {
        ReA.push_back( MatrixOfZeros( NumOfStates, NumOfStates) );
        ImA.push_back( MatrixOfZeros( NumOfStates, NumOfStates) );
    }

    vector<double> MomentaSet;
    for(int q = 0; q < NumOfStates; q++)
    {
        MomentaSet.push_back(  2.*M_PI*(-(NumOfStates - 1)*1./2 + q)/L  );
        cout<< "\n\n" << MomentaSet[q];
    }





}


//R_WF_approx(vector<double> X, double L, double c)
//R_pushed_WF_approx(vector<double> X, double P, double L, double c)
