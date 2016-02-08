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


private:
    vector<vector<int> > Perms;
    int _N = 0;
    vector<int> _PermPool;
    vector<int> _PermSigns;
    bool _Permity = true;


public:


    void SaveVecToFile(vector<double> from,
                       std::ofstream &to);

    void SaveVecFromFile(int veclength,
                         std::ifstream &from,
                         vector<double> &to);

    void SaveMatrixToFile(vector<vector<double>> from,
                          std::ofstream &to);

    void SaveMatrixFromFile(int rows,
                            int columns,
                            std::ifstream &from,
                            vector<vector<double>> &to);




    void setPermSize(int N);
    int  getPermSize();
    void PermSwap(int a, int b);
    void Permute(int k, int size);
    void appendPerm();

    std::complex<double> sProduct(int n, int PermNumber,
                                  double c,
                                  vector<double> X,
                                  vector<double> K);

    std::complex<double> sProduct_unnormalized(int n,
                                              int PermNumber,
                                              double c,
                                              vector<double> X,
                                              vector<double> K);

    std::complex<double> sProduct_strong(int n,
                                         int PermNumber,
                                         double c,
                                         vector<double> X,
                                         vector<double> K);

    int SignFunction(double x1, double x2);

    int GetPermEl(int i, int j);

    std::complex<double> WaveFunction(int n, double c,
                                      vector<double> X,
                                      vector<double> K,
                                      double norm);

    std::complex<double> WaveFunction_unnormalized(int n, double c,
                                      vector<double> X,
                                      vector<double> K,
                                      double norm);

    std::complex<double> WaveFunction_strong(int n,
                                             double c,
                                             vector<double> X,
                                             vector<double> K,
                                             double norm);

    double GetRandom(double min, double max);

    vector<double> ParticlePositions;

    void Positions(vector<double> X0,
                   vector<double> &X1,
                   double delta,
                   int n,
                   double L);

    void AddCollection(vector<vector<double> > &Chain,
                       vector<double>  Set);

    void Metropolis(int steps,
                    int N,
                    double c,
                    vector<vector<double> > &Chain,
                    vector<double> &WaveFValues,
                    vector<double>  &X0,
                    vector<double>  &X1,
                    vector<double>  K,
                    double delta,
                    double L,
                    double norm);

    void Metropolis_unnormalized(int OutSteps,
                                 int InsideSteps,
                                int N,
                                double c,
                                vector<vector<double> > &Chain,
                                vector<double> &WaveFValues,
                                vector<double> &X1,
                                vector<double> K,
                                double delta,
                                double L,
                                double norm);

    void Metropolis_strong(int OutSteps,
                           int InsideSteps,
                           int N,
                           double c,
                           vector<vector<double> > &Chain,
                           vector<double> &WaveFValues,
                           vector<double> &X1,
                           vector<double> K,
                           double delta,
                           double L,
                           double norm);

    int  MaximumPosition(vector<double> WaveFValues);

    vector<double> PhaseCalculation(std::complex<double> WF);

    void MinimaFinder(vector<double> &Minima,
                      vector<vector<double>> &PhaseMatrix,
                      vector<vector<double>> &WFMatrix,
                      vector<vector<double> > &SinMatrix,
                      vector<vector<double> > &CosMatrix,
                      int NumOfBisections,
                      int InitialDivision,
                      int n,
                      int ChainElement,
                      double JumpRestriction,
                      double c,
                      double L,
                      vector<vector<double> > Chain,
                      vector<double> K,
                      double norm);

    void MinimaFinder_unnormalized(vector<double> &Minima,
                                  vector<vector<double>> &PhaseMatrix,
                                  vector<vector<double>> &WFMatrix,
                                  vector<vector<double> > &SinMatrix,
                                  vector<vector<double> > &CosMatrix,
                                  int NumOfBisections,
                                  int InitialDivision,
                                  int n,
                                  int ChainElement,
                                  double JumpRestriction,
                                  double c,
                                  double L,
                                  vector<vector<double> > Chain,
                                  vector<double> K,
                                  double norm);

    void MinimaFinder_strong(vector<double> &Minima,
                                  vector<vector<double>> &PhaseMatrix,
                                  vector<vector<double>> &WFMatrix,
                                  vector<vector<double> > &SinMatrix,
                                  vector<vector<double> > &CosMatrix,
                                  int NumOfBisections,
                                  int InitialDivision,
                                  int n,
                                  int ChainElement,
                                  double JumpRestriction,
                                  double c,
                                  double L,
                                  vector<vector<double> > Chain,
                                  vector<double> K,
                                  double norm);



    void IntegrationCollections(vector<double> K1,
                                vector<double> K2,
                                vector<vector<double>> &AllX,
                                vector<double> &AllWs,
                                vector<double> &AllNorms,
                                double L,
                                int PolynomialDegree,
                                int IntervalDivision);

    void ScalarProduct(vector<double> K1,
                       vector<double> K2,
                       vector<vector<double>> AllX,
                       vector<double> AllWs,
                       vector<double> AllNorms,
                       double c,
                       vector<double> Measured_X,
                       vector<double> &ScalarProducts_re,
                       vector<double> &ScalarProducts_im,
                       double norm1,
                       double norm2);

    void Scalar_GaussLegendre_Loop(int NumOfQuasimomentaSets,
                                               int PolynomialDegree,
                                               int IntervalDivision,
                                               vector<double> Measured_X,
                                               vector<vector<double>> K_Matrix,
                                               vector<double> K2,
                                               vector<double> Norms,
                                               double L,
                                               double c,
                                               double norm,
                                               std::ofstream &sc_re,
                                               std::ofstream &sc_im);

    void Scalar_MonteCarlo_Loop(int Type,
                                int NumOfPositions,
                                int NumOfParticles,
                                int NumOfQuasimomentaSets,
                                vector<double> Measured_X,
                                vector<vector<double>> K_Matrix,
                                vector<double> K2,
                                vector<double> Norms,
                                double L,
                                double c,
                                double norm,
                                std::ofstream &sc_re,
                                std::ofstream &sc_im);


    void RandPos_MonteCarloIntegration(int L,
                                       int NumOfParticles,
                                       int NumOfPositions,
                                       vector<vector<double>> &RandPos_Matrix);

    void MonteCarlo_ScalarProduct(int Type,
                                  int L,
                                  vector<double> K1,
                                  vector<double> K2,
                                  vector<vector<double>> RandPos_Matrix,
                                  double c,
                                  vector<double> Measured_X,
                                  vector<double> &ScalarProducts_re,
                                  vector<double> &ScalarProducts_im,
                                  double norm1,
                                  double norm2);

    std::complex<double> WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix,
                                                      vector<double> ScalarProducts_re,
                                                      vector<double> ScalarProducts_im,
                                                      vector<double> Remaining_X,
                                                      double c,
                                                      double t,
                                                      vector<double> Norms);

    void RelevantCoefficients(vector<double> &ScalarProducts_re,
                              vector<double> &ScalarProducts_im,
                              vector<vector<double>> &K_Matrix,
                              double CutOff_Coefficient,
                              vector<double> &Norms);

    void Metropolis_TimeDependent(int steps,
                                  int N,
                                  double c,
                                  double t,
                                  vector<vector<double> > &Chain,
                                  vector<double> &WaveFValues,
                                  vector<double>  &X0,
                                  vector<double>  &X1,
                                  vector<vector<double>>  K_Matrix,
                                  vector<double> ScalarProducts_re,
                                  vector<double> ScalarProducts_im,
                                  double delta,
                                  double L,
                                  vector<double> Norms);



// ########################## REPULSIVE CASE #########################################################################################################################################################################################################################################################################################################################################

    std::complex<double> R_sProduct(int n,
                                    int PermNumber,
                                    double c,
                                    vector<double> X,
                                    vector<double> K_re,
                                    vector<double> K_im);

    std::complex<double> R_WaveFunction(int n,
                                        double c,
                                        vector<double> X,
                                        vector<double> K_re,
                                        vector<double> K_im,
                                        double normalizacja);


    double R_TwoBodyWF_approx_merged(double x1, double x2,double L, double c);

    double R_WF_approx(vector<double> X, double L, double c);

    std::complex<double> R_pushed_WF_approx(vector<double> X, double P, double L, double c);

    void R_Metropolis_approx(int OutSteps,
                                       int InsideSteps,
                                       int N,
                                       double c,
                                       vector<vector<double> > &Chain,
                                       vector<double> &WaveFValues,
                                       vector<double> &X1,
                                       double delta,
                                       double L
                                       );

    void R_IntegrationCollections(vector<double> K1_re,
                                  vector<double> K2_re,
                                  vector<vector<double>> &AllX,
                                  vector<double> &AllWs,
                                  vector<double> &AllNorms,
                                  double L,
                                  int PolynomialDegree,
                                  int IntervalDivision);

    void R_IntegrationCollections_approx(int N1,
                                         int N2,
                                         vector<vector<double>> &AllX,
                                         vector<double> &AllWs,
                                         vector<double> &AllNorms,
                                         double L,
                                         int PolynomialDegree,
                                         int IntervalDivision);

    void R_ScalarProduct(vector<double> K1_re,
                         vector<double> K1_im,
                         vector<double> K2_re,
                         vector<double> K2_im,
                         const vector<vector<double>>& AllX,
                         const vector<double>& AllWs,
                         const vector<double>& AllNorms,
                         double c,
                         vector<double> Measured_X,
                         vector<double> &ScalarProducts_re,
                         vector<double> &ScalarProducts_im,
                         double norm1,
                         double norm2);


    void R_ScalarProduct_approx(double P,
                                const vector<vector<double> > &AllX,
                                const vector<double> &AllWs,
                                const vector<double> &AllNorms,
                                double c,
                                double L,
                                vector<double> Measured_X,
                                vector<double> &ScalarProducts_re,
                                vector<double> &ScalarProducts_im);

    void R_MC_ScalarProduct_approx( int L,
                                    double c,
                                    double P,
                                    int NumOfParticlesAfterMeasurement,
                                    int NumberOfPoints,
                                    double range,
                                    vector<double> Measured_X,
                                    vector<double> &ScalarProducts_re,
                                    vector<double> &ScalarProducts_im,
                                    vector<double> &NormsAfterMeasurement
                                    );



    void R_RelevantCoefficients(vector<double> &ScalarProducts_re,
                                vector<double> &ScalarProducts_im,
                                vector<vector<double>> &K_Matrix_re,
                                vector<vector<double>> &K_Matrix_im,
                                double CutOff_Coefficient);

    void R_MonteCarlo_ScalarProduct(int Type,
                                    int L,
                                    vector<double> K1_re,
                                    vector<double> K1_im,
                                    vector<double> K2_re,
                                    vector<double> K2_im,
                                    const vector<vector<double> > &RandPos_Matrix,
                                    double c,
                                    vector<double> Measured_X,
                                    vector<double> &ScalarProducts_re,
                                    vector<double> &ScalarProducts_im,
                                    double norm1,
                                    double norm2);

    std::complex<double> R_WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix_re,
                                                        vector<vector<double> > K_Matrix_im,
                                                        vector<double> ScalarProducts_re,
                                                        vector<double> ScalarProducts_im,
                                                        vector<double> Remaining_X,
                                                        double c,
                                                        double phi,
                                                        double t,
                                                        vector<double> Norms);

    void R_Metropolis_TimeDependent(int N,
                                    double c,
                                    double t,
                                    vector<vector<double> > &Chain,
                                    vector<double> &WaveFValues,
                                    vector<double>  &X0,
                                    vector<double>  &X1,
                                    vector<vector<double>>  K_Matrix_re,
                                    vector<vector<double>>  K_Matrix_im,
                                    vector<double> ScalarProducts_re,
                                    vector<double> ScalarProducts_im,
                                    double delta,
                                    double L,
                                    double phi,
                                    vector<double> Norms);



//signals:

//public slots:
};


#endif // PERMUTATOR_H
