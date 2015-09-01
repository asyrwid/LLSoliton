#ifndef PERMUTATOR_H
#define PERMUTATOR_H
#include <complex>
#include <vector>
#include <iostream>
#include <omp.h>
using std::vector;
using std::cout;

class permutator
{

public:
    // explicit Permutator();
//    ~Permutator();


private:
    vector<vector<int> > Perms;
    int _N = 0;
    vector<int> _PermPool;
    vector<int> _PermSigns;
    bool _Permity = true;


public:
    void setPermSize(int N);
    int  getPermSize();
    void PermSwap(int a, int b);
    void Permute(int k, int size);
    void appendPerm();

    std::complex<double> sProduct(int n, int PermNumber, double c, vector<double> X, vector<double> K);
    int SignFunction(double x1, double x2);
    int GetPermEl(int i, int j);
    std::complex<double> WaveFunction(int n, double c, vector<double> X, vector<double> K, double norm);
    double GetRandom(double min, double max);
    vector<double> ParticlePositions;
    void Positions(vector<double> X0, vector<double> &X1, double delta, int n, double L);
    void AddCollection(vector<vector<double> > &Chain, vector<double>  Set);

    void Metropolis(int N,double c, vector<vector<double> > &Chain, vector<double> &WaveFValues,  vector<double>  &X0,vector<double>  &X1, vector<double>  K, double delta, double L, double norm);
    int  MaximumPosition(vector<double> WaveFValues);
    vector<double> PhaseCalculation(std::complex<double> WF);
    void MinimaFinder(vector<double> &Minima, vector<vector<double>> &PhaseMatrix, vector<vector<double>> &WFMatrix, vector<vector<double> > &SinMatrix, vector<vector<double> > &CosMatrix, int NumOfBisections, int InitialDivision, int n, int ChainElement, double JumpRestriction ,double c, double L, vector<vector<double> > Chain, vector<double> K, double norm);

    void IntegrationCollections(vector<double> K1, vector<double> K2, vector<vector<double>> &AllX, vector<double> &AllWs, vector<double> &AllNorms, double L, int PolynomialDegree, int IntervalDivision);
    void ScalarProduct(vector<double> K1, vector<double> K2, vector<vector<double>> AllX, vector<double> AllWs, vector<double> AllNorms, double c, vector<double> Measured_X, vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im, double norm1, double norm2);
    std::complex<double> WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, vector<double> Remaining_X, double c, double t, vector<double> Norms);
    void RelevantCoefficients(vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im, vector<vector<double>> &K_Matrix, double CutOff_Coefficient, vector<double> &Norms);
    void Metropolis_TimeDependent(int N,double c, double t, vector<vector<double> > &Chain, vector<double> &WaveFValues,  vector<double>  &X0,vector<double>  &X1, vector<vector<double>>  K_Matrix, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, double delta, double L, vector<double> Norms);



// ########################## REPULSIVE CASE #########################################################################################################################################################################################################################################################################################################################################

    std::complex<double> R_sProduct(int n, int PermNumber, double c, vector<double> X, vector<double> K_re, vector<double> K_im);
    std::complex<double> R_WaveFunction(int n, double c, vector<double> X, vector<double> K_re, vector<double> K_im);

    void R_IntegrationCollections(vector<double> K1_re, vector<double> K2_re, vector<vector<double>> &AllX, vector<double> &AllWs, vector<double> &AllNorms, double L, int PolynomialDegree, int IntervalDivision);
    void R_ScalarProduct(vector<double> K1_re, vector<double> K1_im, vector<double> K2_re, vector<double> K2_im, vector<vector<double>> AllX, vector<double> AllWs, vector<double> AllNorms, double c, vector<double> Measured_X, vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im);
    void R_RelevantCoefficients(vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im, vector<vector<double>> &K_Matrix_re, vector<vector<double>> &K_Matrix_im, double CutOff_Coefficient);
    std::complex<double> R_WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix_re, vector<vector<double> > K_Matrix_im, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, vector<double> Remaining_X, double c, double t);
    void R_Metropolis_TimeDependent(int N,double c, double t, vector<vector<double> > &Chain, vector<double> &WaveFValues,  vector<double>  &X0,vector<double>  &X1, vector<vector<double>>  K_Matrix_re, vector<vector<double>>  K_Matrix_im, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, double delta, double L);



//signals:

//public slots:
};
#endif // PERMUTATOR_H
