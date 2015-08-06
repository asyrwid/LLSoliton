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
    std::complex<double> WaveFunction(int n, double c, vector<double> X, vector<double> K);
    double GetRandom(double min, double max);
    vector<double> ParticlePositions;
    void Positions(vector<double> X0, vector<double> &X1, double delta, int n, double L);
    void AddCollection(vector<vector<double> > &Chain, vector<double>  Set);
    void Metropolis(int N,double c, vector<vector<double> > &Chain,  vector<double> &WaveFValues, vector<double>  &X0, vector<double>  &X1, vector<double>  K, double delta, int n, double L);
    int  MaximumPosition(vector<double> WaveFValues);
    vector<double> PhaseCalculation(std::complex<double> WF);
    void MinimaFinder(vector<double> &Minima, vector<vector<double>> &PhaseMatrix, vector<vector<double>> &WFMatrix, int NumOfBisections, int InitialDivision, int n, int ChainElement, double JumpRestriction, double c, double L, vector<vector<double> > Chain, vector<double> K);

//signals:

//public slots:
};
#endif // PERMUTATOR_H
