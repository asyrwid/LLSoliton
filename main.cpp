#include "permutator.h"
#include "decomposition.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <omp.h>        /* OpenMP */
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <list>



typedef std::complex<double> Complex;

# define M_PI  3.14159265358979323846  /* pi definition */
int main()
{
    double t0=clock();
    srand( time( NULL ) );


//===========================================================================================================
//*************** INITIAL VALUES AND OBJECTS ****************************************************************
//===========================================================================================================

    int    Steps = 2000;  // number of metropolis steps
    int    N = 6;         // number of particles
    double c = 1.;      // coupling constant
    double L = 1.0;       // system length
    double delta = 0.2;   // maximal "jump" between single particle position from one realization to another
    double norm = 12.1486;// norm of initial state

    // initial collection of particle positions
    vector<double> X0 = {.1,.25,.3,.5,.8,.9};
    vector<double> X1 = X0;

    // quasimomentas
    vector<double> K ={3.90183, 5.49048, 7.05365, 11.7959, 13.3591, 14.9477};

    // chain of position
    vector<vector<double> > Chain;

    //chain of wave function values
    vector<double> WaveFValues;

    //vector of wave function minima and Phase/WF (+ Sin and Cos of Phase) PlotPoints
    vector<double> Minima;
    vector<vector<double> > PhaseMatrix;
    vector<vector<double> > SinMatrix;
    vector<vector<double> > CosMatrix;
    vector<vector<double> > WFMatrix;

permutator _Permutator;
decomposition _Decomposition;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //cout << _Permutator.getPermSize(); // initially _N=0;
    _Permutator.setPermSize(N); // de facto particle number
    //cout << _Permutator.getPermSize(); // N is established now!
    _Permutator.Permute(N,N); // all permutations preparation
    _Permutator.AddCollection(Chain, X0);

    // Num of quasimomenta after measurement + num of quasimomenta collections
    int NumOfQuasimomenta = 3;
    int NumOfQuasimomentaSets = 165;

    std::ifstream norms;
    norms.open("E:/Projekty cpp/Permutator/Normy_3_c10.txt");
    vector<double> Norms;
    _Permutator.SaveVecFromFile( NumOfQuasimomentaSets, norms, Norms);
    norms.close();

    _Decomposition.Coefficients_CramerMethod(5,{.5},17,L);

// ================ SCALAR PRODUCTS CALCULATION and saving to files =============================================================
/*


    std::ifstream k_matrix ("E:/Projekty cpp/Permutator/Kwazi_3_c10.txt");
    vector<vector<double>> K_Matrix;

    _Permutator.SaveMatrixFromFile(NumOfQuasimomentaSets,NumOfQuasimomenta,k_matrix,K_Matrix);
    k_matrix.close();


    std::ofstream chain_strong ("E:/projekty cpp/Permutator/chain_s.txt");
    std::ofstream phaseMat ("E:/projekty cpp/Permutator/phaseMat.txt");
    std::ofstream wfMat  ("E:/projekty cpp/Permutator/wfMat.txt");
    std::ofstream sinMat ("E:/projekty cpp/Permutator/sinMat.txt");
    std::ofstream cosMat ("E:/projekty cpp/Permutator/cosMat.txt");


    int OutSteps = 1;
    int InsideSteps = 1000;

    _Permutator.Metropolis_unnormalized(OutSteps, InsideSteps, N, c, Chain, WaveFValues, X1, K, delta, L, 1.);

    cout<< "\n\n" << Chain.size() << "\n";

    int NumOfBisections = 10;
    int InitialDivision = 100;
    double JumpRestriction = .1;

    for( int ChainElement = 0; ChainElement< Chain.size(); ChainElement++)
    {
        _Permutator.MinimaFinder_unnormalized(Minima, PhaseMatrix, WFMatrix, SinMatrix, CosMatrix,
                                    NumOfBisections, InitialDivision, N, ChainElement,
                                    JumpRestriction, c, L, Chain, K, 1.);

        if(ChainElement%100 == 0)
        {
           cout<< "\r"<< ChainElement << " ";
        }
    }


    _Permutator.SaveMatrixToFile(Chain, chain_strong);
    _Permutator.SaveMatrixToFile(PhaseMatrix, phaseMat);
    _Permutator.SaveMatrixToFile(WFMatrix, wfMat);
    _Permutator.SaveMatrixToFile(SinMatrix, sinMat);
    _Permutator.SaveMatrixToFile(CosMatrix, cosMat);





/*
    int NumOfPoints = 100;
    std::ofstream wfs ("E:/projekty cpp/Permutator/wfki.txt");
    std::ofstream args ("E:/projekty cpp/Permutator/argsy.txt");

    for(int j = 0; j < NumOfPoints; j++)
    {
        vector<double> X = X0;
        X.push_back(1./NumOfPoints * j);

        double hr =pow( fabs(_Permutator.WaveFunction_strong(N, c, X, K, 1.)),2.);
        double argi = arg (_Permutator.WaveFunction_strong(N, c, X, K, 1.));

        wfs << hr << " ";
        args << argi << " ";
    }
*/

/*
    std::ofstream chain_approx ("E:/projekty cpp/Permutator/chain_app.txt");

    int OutSteps = 1;
    int InsideSteps = 5000000;


    _Permutator.R_Metropolis_approx(OutSteps, InsideSteps, N, c, Chain, WaveFValues, X1, delta, L );

    _Permutator.SaveMatrixToFile(Chain, chain_approx);





// ================= SCALAR PRODUCT (MONTE CARLO) ===============================================================================
/*
    vector<vector<double>> RandPos_Matrix;
    int NumOfParticles = 3;
    int NumOfPositions = 500000;

    _Permutator.RandPos_MonteCarloIntegration(L, NumOfParticles, NumOfPositions, RandPos_Matrix);

    std::ofstream sc_re ("E:/projekty cpp/Permutator/Scalar_re_3_mc.txt");
    std::ofstream sc_im ("E:/projekty cpp/Permutator/Scalar_im_3_mc.txt");
    vector<double> K2 = K;
    vector<double> Measured_X ={.1,.3,.5};

    _Permutator.Scalar_MonteCarlo_Loop(2, NumOfPositions, NumOfParticles, NumOfQuasimomentaSets,
                           Measured_X, K_Matrix, K2, Norms, L, c, norm, sc_re, sc_im);

*/
//============== SCALAR PRODUCT (GAUSS-LEGENDRE)  ===============================================================================
/*
    int PolynomialDegree = 5;
    int IntervalDivision = 12;

    vector<double> K2 = K;
    vector<double> Measured_X = {.1,.3,.5};

    std::ofstream sc_re ("E:/projekty cpp/Permutator/Scalar_re_3_c10.txt");
    std::ofstream sc_im ("E:/projekty cpp/Permutator/Scalar_im_3_c10.txt");

    _Permutator.Scalar_GaussLegendre_Loop(NumOfQuasimomentaSets,
                                          PolynomialDegree,
                                          IntervalDivision,
                                          Measured_X, K_Matrix, K2,
                                          Norms, L, c, norm, sc_re, sc_im);



// ================================ TIME EVOLUTION ================================================================
/*
    //std::ifstream k_matrix  ("E:/Projekty cpp/Permutator/Kwazi_3_c1.txt");
    std::ifstream scprod_re ("E:/projekty cpp/Permutator/Scalar_re_3_c10.txt");
    std::ifstream scprod_im ("E:/projekty cpp/Permutator/Scalar_im_3_c10.txt");

    //vector<vector<double>> K_Matrix;
    vector<double> ScalarProducts_re, ScalarProducts_im;

    //_Permutator.SaveMatrixFromFile(NumOfQuasimomentaSets,NumOfQuasimomenta,k_matrix,K_Matrix);
    _Permutator.SaveVecFromFile(NumOfQuasimomentaSets,scprod_re,ScalarProducts_re);
    _Permutator.SaveVecFromFile(NumOfQuasimomentaSets,scprod_im,ScalarProducts_im);
    //k_matrix.close();
    scprod_re.close();
    scprod_im.close();

// ======== WE CHOOSE RELEVANT ELEMENTS ===========================================================================

/*
    double CutOff_Coefficient =  0.0005;

    _Permutator.RelevantCoefficients(ScalarProducts_re, ScalarProducts_im, K_Matrix, CutOff_Coefficient, Norms);

    cout << "Number of relevant elements: " << ScalarProducts_im.size() << "\n";

// ========= METROPOLIS FOR DIFFERENT TIME MOMENTS ================================================================


    int Steps_t = 5000;
    int N_Remaining = (K_Matrix.at(0)).size();

    vector<vector<vector<double>>> Chain_T_dependence;

    std::ofstream chain_t ("E:/projekty cpp/Permutator/X_T.txt");

    vector<vector<double>> Chain_t;

    for(int j = 0; j < 250; j++)
    {
        double t = 0.02 * j;

        cout << "\n\n" << "time value: " << t << "\n\n";

        //vector<vector<double>> Chain_t;
        vector<double> WaveFValues_t;
        vector<double> X0_t;
        for(int i = 0; i < N_Remaining; i++)
        {
            X0_t.push_back(_Permutator.GetRandom(0.,1.));
        }
        vector<double> X1_t = X0_t;
        _Permutator.AddCollection(Chain_t, X0_t);


        _Permutator.Metropolis_TimeDependent(Steps_t, N_Remaining, c, t, Chain_t,
                                             WaveFValues_t, X0_t, X1_t, K_Matrix,
                                             ScalarProducts_re, ScalarProducts_im, delta,
                                             L, Norms);

        for(unsigned s = 0; s < Chain_t.size(); s++)
        {
            for(int y = 0; y < N_Remaining; y++ )
            {
                chain_t << ( Chain_t.at(s) ).at(y) << " ";
            }
            chain_t << "\n";
        }
        chain_t << "\n";
        Chain_T_dependence.push_back(Chain_t);
        Chain_t.clear();
    }
    chain_t.close();





/* ################################################################################################################
**************************************** NORMALIZATION CHECK ******************************************************
################################################################################################################ */


/*
    int PolynomialDegree = 5;
    int IntervalDivision = 11;

    vector<vector<double>> AllX;
    vector<double> AllWs, AllNorms;
    vector<double> ScalarProducts_re, ScalarProducts_im, NormsAfterMeasurement;


    int NumOfParticlesAfterMeasurement = 10;
    int NumberOfPoints = 10000000;
    double  range = 1./fabs(c) *40.;



    vector<double> Measured_X={.4};

    //_Permutator.R_IntegrationCollections_approx(4,4, AllX, AllWs, AllNorms,
    //                                            L, PolynomialDegree, IntervalDivision);

    for(int j = 0; j < 5; j++)
    {
        double P = 2. * (j-2.) * M_PI / L;


        //_Permutator.R_ScalarProduct_approx(P, AllX, AllWs, AllNorms, c, L,
        //                                   Measured_X, ScalarProducts_re, ScalarProducts_im);


        _Permutator.R_MC_ScalarProduct_approx(L, c, P, NumOfParticlesAfterMeasurement,
                                              NumberOfPoints, range, Measured_X,
                                              ScalarProducts_re, ScalarProducts_im, NormsAfterMeasurement);

    }

    std::ofstream scre ("C:/Users/Andrzej Syrwid/Desktop/scre.txt");
    std::ofstream scim ("C:/Users/Andrzej Syrwid/Desktop/scim.txt");

    _Permutator.SaveVecToFile(ScalarProducts_re,scre);
    _Permutator.SaveVecToFile(ScalarProducts_im,scim);


/*


    int NumOfPoints = 100;
    double WF6, WF3;
    double Wf6Value = 0.;
    double Wf3Value = 0.;

    std::ofstream w6 ("E:/projekty cpp/Permutator/wf6.txt");
    std::ofstream w3 ("E:/projekty cpp/Permutator/wf3.txt");


    for(int j = 0; j < NumOfPoints ; j++)
    {
        vector<double> X_01 = {.4,.39, .41};
        vector<double> X_02 = { .39, .41};

        X_01.push_back( 1./NumOfPoints * j );
        X_02.push_back( 1./NumOfPoints * j );


        Complex wf_ped(0., 0.);
        for(int s = 0; s < 11; s++)
        {
            double p = 2. * (j-5.) * M_PI / L;

            wf_ped = wf_ped + _Permutator.R_pushed_WF_approx(X_02, p, L, c) * Complex(ScalarProducts_re[s], ScalarProducts_im[s]) ;
        }


        WF6 = pow( fabs ( _Permutator.R_WF_approx(X_01, L, c) ), 2.);
        WF3 = pow( fabs ( wf_ped ), 2.);

        Wf6Value = Wf6Value + WF6 * 1./NumOfPoints;
        Wf3Value = Wf3Value + WF3 * 1./NumOfPoints;

        cout << Wf6Value << " , " << Wf3Value << " , " << fabs(Wf6Value - Wf3Value) << "\n";

        w6 << WF6 << " ";
        w3 << WF3 << " ";

        X_01.clear();
        X_02.clear();
    }



  /*
    int NumOfPoints = 100;
    double WF6, WF3;
    double Wf6Value = 0.;
    double Wf3Value = 0.;

    std::ofstream w6 ("E:/projekty cpp/Permutator/wf6.txt");
    std::ofstream w3 ("E:/projekty cpp/Permutator/wf3.txt");

    for(int j = 0; j < NumOfPoints ; j++)
    {
        vector<double> X_01 = {.1,.3,.5,.9,.22};
        vector<double> X_02 = { .9,.22};

        X_01.push_back( 1./NumOfPoints * j );
        X_02.push_back( 1./NumOfPoints * j );

        WF6 = pow( fabs ( _Permutator.WaveFunction(N, c, X_01, K, norm) ), 2.);
        WF3 = pow( fabs ( _Permutator.WaveFunctionAfterMeasurement(K_Matrix, ScalarProducts_re, ScalarProducts_im, X_02,  c, 0., Norms) ), 2.);

        Wf6Value = Wf6Value + WF6 * 1./NumOfPoints;
        Wf3Value = Wf3Value + WF3 * 1./NumOfPoints;

        cout << Wf6Value << " , " << Wf3Value << " , " << fabs(Wf6Value - Wf3Value) << "\n";

        w6 << WF6 << " ";
        w3 << WF3 << " ";

        X_01.clear();
        X_02.clear();
    }





*/












/* ================================================================================================================
###################################################################################################################
----------------------------- REPULSIVE CASE ----------------------------------------------------------------------
###################################################################################################################
================================================================================================================ */

// ===================================================================================================
// ============================= INITIAL STATE =======================================================
// ===================================================================================================

/*

    vector<double> K_re = {9.00091989557964766486923818139542774474648 * pow(10.,-9.),4.2364030224662437344956920221017426797516 * pow(10.,-9.), -6.5441527366687672705831361336104023962434 * pow(10.,-9.),-6.69317018137712412878179406988676802825466 * pow(10.,-9.)};
    vector<double> K_im = {-0.77847205566101472282729215878601769487839558505537 , -2.4441586119425518633592182578634431956670400966531, \
                           2.4441586117379878391700426542615305397009649026454, 0.77847205586557874701646776238793035084447077906314};


/*
    vector<double> R_X = {.1, .2, .223, .14};

    std::complex<double> R_Wave;
    R_Wave = _Permutator.R_WaveFunction(N , c, R_X, K_re, K_im, 685.715536036939239371756525487573480235208398835);
    cout<< R_Wave.real() << ", " << R_Wave.imag();
*/

// ===================================================================================================
// ================== FILES READING (quasimomenta matrices and norms) ================================
// ===================================================================================================
/*
    std::ifstream k_matrix_re;
    std::ifstream k_matrix_im;
    std::ifstream norms_tc;
    k_matrix_re.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\K_re_3_TC.txt");
    k_matrix_im.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\K_im_3_TC.txt");
    norms_tc.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\Normy_3_TC.txt");

    vector<double> Norms_TC;

    int NumOfQuasimomenta = 3;
    int NumOfQuasimomentaSets = 41;

    vector<vector<double>> K_Matrix_re, K_Matrix_im;
    vector<double> TmpVec_re, TmpVec_im;
    double tmp_re, tmp_im, y;

    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        for(int i = 0; i < NumOfQuasimomenta; i++)
        {
             k_matrix_re >> tmp_re;
             k_matrix_im >> tmp_im;
             TmpVec_re.push_back(tmp_re);
             TmpVec_im.push_back(tmp_im);
        }
        K_Matrix_re.push_back( TmpVec_re );
        K_Matrix_im.push_back( TmpVec_im );
        TmpVec_re.clear();
        TmpVec_im.clear();

        norms_tc >> y;
        Norms_TC.push_back(y);
        //cout << Norms_TC.back() << " " ;
    }

    k_matrix_re.close();
    k_matrix_im.close();
    norms_tc.close();

/*
    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        for(int i = 0; i < NumOfQuasimomenta; i++)
        { cout << K_Matrix_re[j][i] << ", "; }
        cout << "\n";
    }
*/


//====================================================================================================================================
//================= Scalar products calculation (Gauss-Legendre)  ====================================================================
//====================================================================================================================================


/*
    vector<vector<double>> AllX;
    vector<double> AllWs, AllNorms;

    int PolynomialDegree = 5;
    int IntervalDivision = 10;

    vector<double> Measured_X = {.5};

    //vector<double> K_re = K_Matrix_re.at(15);
    //vector<double> K_im = K_Matrix_im.at(15);

    _Permutator.R_IntegrationCollections(K_Matrix_re[0], K_re ,AllX ,AllWs ,AllNorms ,L, PolynomialDegree, IntervalDivision);


    std::ofstream scalar_re ("Scalar_re_3_TC.txt");
    std::ofstream scalar_im ("Scalar_im_3_TC.txt");
    std::fstream sc_re;
    std::fstream sc_im;
    sc_re.open ("Scalar_re_3_TC.txt");
    sc_im.open ("Scalar_im_3_TC.txt");

    vector<double>  ScalarProducts_re, ScalarProducts_im, Sc_re_tmp, Sc_im_tmp;
    vector<int> Parallel_Order;
    Sc_re_tmp = vector<double> (NumOfQuasimomentaSets);
    Sc_im_tmp = vector<double> (NumOfQuasimomentaSets);


    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Parallel for loop
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int u = 0;
    int Par_i1; // iterator in the first parallelized loop
    double Par_norm1; // norm of tha Par_i1-th state of states considered after measurement
    double Par_norm2 = 44.5827500877363236765351775573911399118111440835;
    vector<double> K1_re;
    vector<double> K1_im;

    #pragma omp parallel private(Par_i1 ,K1_re, K1_im, Par_norm1 )
    #pragma omp for
    for(Par_i1 = 0; Par_i1 < NumOfQuasimomentaSets; Par_i1++)
    {
        K1_re = K_Matrix_re.at(Par_i1);
        K1_im = K_Matrix_im.at(Par_i1);
        Par_norm1 = Norms_TC.at(Par_i1);

        _Permutator.R_ScalarProduct(K1_re, K1_im, K_re, K_im, AllX, AllWs, AllNorms, c,
                                    Measured_X, ScalarProducts_re,ScalarProducts_im, Par_norm1, Par_norm2);

        Parallel_Order.push_back( Par_i1 );
        u++;
        if( u % 1 == 0)
        {
            cout<< "\n\n" << u << " "  << " " <<  "( " << ScalarProducts_re.back() << " , " << ScalarProducts_im.back() << " )" << ",     " << pow( ScalarProducts_re.back(), 2.)+pow( ScalarProducts_im.back(), 2.)  << "\n\n" ;;

        }

        //cout << "\n \n Row number of K_Matrix: " << Par_i1 << ", value of re part: " << ScalarProducts_re.back() << ", par order: " << Parallel_Order.back() << "\n \n";
    }
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Sort and writing to file :)
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        int position = Parallel_Order.at(j);

        Sc_re_tmp[ position ] =  ScalarProducts_re.at( j );
        Sc_im_tmp[ position ] =  ScalarProducts_im.at( j );
    }
    ScalarProducts_re = Sc_re_tmp;
    ScalarProducts_im = Sc_im_tmp;

    for(int s = 0; s < NumOfQuasimomentaSets; s++)
    {
        sc_re << ScalarProducts_re.at(s) << " ";
        sc_im << ScalarProducts_im.at(s) << " ";
    }
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    sc_re.close();
    sc_im.close();


/*

    vector<double> K1_re, K2_re;
    vector<double> K1_im, K2_im;
    double Par_norm1; // norm of tha Par_i1-th state of states considered after measurement
    double Par_norm2 = 4.00127409578059014134192414579916474735379 * pow(10., 15.);


    K1_re = K_Matrix_re.at(20);
    K1_im = K_Matrix_im.at(20);
    Par_norm1 = Norms_TC.at(20);
    //K2_re = K_Matrix_re.at(20);
    //K2_im = K_Matrix_im.at(20);
    //Par_norm2 = Norms_TC.at(20);

    _Permutator.R_ScalarProduct(K1_re, K1_im, K_re, K_im, AllX, AllWs, AllNorms, c, Measured_X, ScalarProducts_re,ScalarProducts_im, Par_norm1,Par_norm2);
     cout << " " <<  "( " << ScalarProducts_re.back() << " , " << ScalarProducts_im.back() << " )" << "\n\n" ;




/*====================================================================================================================================
=================== Scalar products calculation (Monte Carlo) ========================================================================
====================================================================================================================================*/

//========================== INITIAL VALUES/ELEMENTS =================================================================================

/*

    vector<vector<double>> RandPos_Matrix;

    int NumOfParticles = 3;
    int NumOfPositions = 10000000;

    _Permutator.RandPos_MonteCarloIntegration(L, NumOfParticles, NumOfPositions, RandPos_Matrix);

    std::ofstream scalar_re ("Scalar_re_3_TC.txt");
    std::ofstream scalar_im ("Scalar_im_3_TC.txt");
    std::fstream sc_re;
    std::fstream sc_im;
    sc_re.open ("Scalar_re_3_TC.txt");
    sc_im.open ("Scalar_im_3_TC.txt");


    vector<double> Measured_X = {.5};

    vector<double>  ScalarProducts_re, ScalarProducts_im, Sc_re_tmp, Sc_im_tmp;
    vector<int> Parallel_Order;
    Sc_re_tmp = vector<double> (NumOfQuasimomentaSets);
    Sc_im_tmp = vector<double> (NumOfQuasimomentaSets);
//====================================================================================================================================

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Parallel for loop
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int u = 0;
    int Par_i1; // iterator in the first parallelized loop
    double Par_norm1; // norm of tha Par_i1-th state of states considered after measurement
    double Par_norm2 = 4.00127409578059014134192414579916474735379 * pow(10., 15.);;
    vector<double> K1_re, K2_re;
    vector<double> K1_im, K2_im;
/*

    K1_re = K_Matrix_re.at(0);
    K1_im = K_Matrix_im.at(0);
    Par_norm1 = Norms_TC.at(0);
    K2_re = K_Matrix_re.at(0);
    K2_im = K_Matrix_im.at(0);
    Par_norm2 = Norms_TC.at(0);

     _Permutator.R_MonteCarlo_ScalarProduct(2,L, K1_re, K1_im, K2_re, K2_im, RandPos_Matrix, -15., Measured_X, ScalarProducts_re, ScalarProducts_im, Par_norm1, Par_norm2);
     cout << " " <<  "( " << ScalarProducts_re.back() << " , " << ScalarProducts_im.back() << " )" << "\n\n" ;




    #pragma omp parallel private(Par_i1 ,K1_re, K1_im, Par_norm1 )
    #pragma omp for
    for(Par_i1 = 0; Par_i1 < NumOfQuasimomentaSets; Par_i1++)
    {
        K1_re = K_Matrix_re.at(Par_i1);
        K1_im = K_Matrix_im.at(Par_i1);
        Par_norm1 = Norms_TC.at(Par_i1);

        _Permutator.R_MonteCarlo_ScalarProduct(1,L, K1_re, K1_im, K_re, K_im, RandPos_Matrix, c, Measured_X, ScalarProducts_re, ScalarProducts_im, Par_norm1, Par_norm2);

        Parallel_Order.push_back( Par_i1 );
        u++;
        if( u % 1 == 0)
        {
            cout << u << " " <<  "( " << ScalarProducts_re.back() << " , " << ScalarProducts_im.back() << " )" << "\n\n" ;
        }

        //cout << "\n \n Row number of K_Matrix: " << Par_i1 << ", value of re part: " << ScalarProducts_re.back() << ", par order: " << Parallel_Order.back() << "\n \n";
    }
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Sort and writing to file :)
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        int position = Parallel_Order.at(j);

        Sc_re_tmp[ position ] = ScalarProducts_re.at( j );
        Sc_im_tmp[ position ] = ScalarProducts_im.at( j );
    }
    ScalarProducts_re = Sc_re_tmp;
    ScalarProducts_im = Sc_im_tmp;

    for(int s = 0; s < NumOfQuasimomentaSets; s++)
    {
        sc_re << ScalarProducts_re.at(s) << " ";
        sc_im << ScalarProducts_im.at(s) << " ";
    }
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    sc_re.close();
    sc_im.close();





/*====================================================================================================================================
====================================================================================================================================*/
/*

    std::ifstream scprod_re;
    std::ifstream scprod_im;

    scprod_re.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\Scalar_re_3_TC.txt");
    scprod_im.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\Scalar_im_3_TC.txt");



    vector<double>  ScalarProducts_re, ScalarProducts_im;
    double tmp_sc_re, tmp_sc_im;

    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        scprod_re >> tmp_sc_re;
        scprod_im >> tmp_sc_im;

        ScalarProducts_re.push_back(tmp_sc_re);
        ScalarProducts_im.push_back(tmp_sc_im);
    }

    scprod_re.close();
    scprod_im.close();




    int NumOfPoints = 200;
    double WF6, WF3;
    double Wf6Value = 0.;
    double Wf3Value = 0.;


    std::ofstream wf6 ("wf6.txt");
    std::ofstream wf3 ("wf3.txt");

    std::fstream w6;
    std::fstream w3;

    w6.open ("wf6.txt");
    w3.open ("wf3.txt");
    double Initialnorm = 44.5827500877363236765351775573911399118111440835;
    //double Initialnorm = Norms_TC.at(11);
    //vector<double> K_re = K_Matrix_re.at(11);
    //vector<double> K_im = K_Matrix_im.at(11);



    for(int j = 0; j < NumOfPoints ; j++)
    {
        vector<double> X_01 = {.5,.51,.49};
        vector<double> X_02 = { .51,.49};

        X_01.push_back( 1./NumOfPoints * j );
        X_02.push_back( 1./NumOfPoints * j );

        WF6 = pow( fabs ( _Permutator.R_WaveFunction(N, c, X_01, K_re, K_im, Initialnorm ) ), 2.);
        WF3 = pow( fabs ( _Permutator.R_WaveFunctionAfterMeasurement(K_Matrix_re, K_Matrix_im,  ScalarProducts_re,  ScalarProducts_im, X_02, c,0.,0., Norms_TC) ), 2.);

        Wf6Value = Wf6Value + WF6 * 1./NumOfPoints;
        Wf3Value = Wf3Value + WF3 * 1./NumOfPoints;

        cout << Wf6Value << " , " << Wf3Value << " , " << fabs(Wf6Value - Wf3Value) << "\n";

        w6 << WF6 << " ";
        w3 << WF3 << " ";

        X_01.clear();
        X_02.clear();

    }




/*====================================================================================================================================
====================================================================================================================================*/

double t1=clock();
cout<< "\n\n\n"<< (t1-t0)/CLOCKS_PER_SEC << "\n\n";

}



