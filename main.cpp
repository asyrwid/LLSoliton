#include "permutator.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <omp.h>        /* OpenMP */
#include <cstdlib>
#include <cstdio>
#include <ctime>




# define M_PI  3.14159265358979323846  /* pi definition */
int main()
{
    double t0=clock();
    srand( time( NULL ) );


//===========================================================================================================
//***********************************************************************************************************
//*************** INITIAL VALUES AND OBJECTS ****************************************************************
//***********************************************************************************************************
//===========================================================================================================

    int    Steps = 2000; // number of metropolis steps
    int    N = 6; // number of particles
    double c = 0.08; // coupling constant
    double L = 1.0; // system length
    double delta = 0.05; // maximal "jump" between single particle position from one realization to another

    // initial collection of particle positions
    vector<double> X0 = {0.5,.1,.15,.2,.25,.3,.35,.4};
    vector<double> X1 = X0;

    // quasimomentas
    vector<double> K ={-0.560631, -0.0750838, 0.410172, 5.87301, 6.35827, 6.84382};

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

//***********************************************************************************************************

permutator _Permutator;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //cout << _Permutator.getPermSize(); // initially _N=0;
    _Permutator.setPermSize(N); // de facto particle number
    //cout << _Permutator.getPermSize(); // N is established now!
    _Permutator.Permute(N,N); // all permutations preparation
    _Permutator.AddCollection(Chain, X0);


    std::ifstream norms;
    norms.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\normy_3_2.txt");

    vector<double> Norms;

    for(int i = 0 ; i < 969; i++)
    {
        double y;
        norms >> y;
        Norms.push_back(y);
        //cout << Norms.back() << " " ;
    }

    double norm = 1.;

    norms.close();

    cout<< "\n\n\n" << fabs( _Permutator.WaveFunction(N, c, X0, K, norm) ) << "\n\n" ;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/*


// ==================== METROPOLIS ==========================================================================
// ==========================================================================================================

    for(int i = 0; i < Steps; i++)
    {
    _Permutator.Metropolis(N,c,Chain,WaveFValues,X0,X1,K,delta,L);
        if (i % 10 == 0)
        {
            cout << i << ", ";
        }
    }

    //for(unsigned i = 0; i < WaveFValues.size(); i++)
    //{
    //cout << WaveFValues.at(i);
    //}

// %%%%%%%%%% writing data from Metropolis %%%%%%%%%%
    std::ofstream create01 ("x.txt");
    std::fstream fs;
    fs.open ("x.txt");

    for(unsigned i = 0; i < Chain.size(); i++)
    {
        for(int j = 0; j < N; j++ )
        {
            fs << (Chain.at(i)).at(j) << " ";
        }
        fs << "\n";
    }
    fs.close();

    std::ofstream create02 ("wf.txt");
    fs.open ("wf.txt");

    for(unsigned i = 0; i < WaveFValues.size(); i++)
    {
            fs << WaveFValues.at(i) << " ";
    }
    fs << "\n";
    fs.close();

// %%%%%%%%%% maximum of probability %%%%%%%%%%
    int MaxPos = _Permutator.MaximumPosition(WaveFValues);
    X0 = Chain.at(MaxPos);
    for(int y = 0; y < N; y++)
    {
        cout << X0[y];
    }




*/



// ****** LAST PROBABILITY AND PHASE OF WAVE FUNCTION *******************************************************
//***********************************************************************************************************
/*
// %%%%%%%%%% txt files preparation %%%%%%%%%%
    std::ofstream create1 ("WFS.txt");
    std::ofstream create2 ("WFPhase.txt");
    std::ofstream create3 ("WFSin.txt");
    std::ofstream create4 ("WFCos.txt");
    std::fstream fs1;
    std::fstream fs2;
    std::fstream fs3;
    std::fstream fs4;
    fs1.open ("WFS.txt");
    fs2.open ("WFPhase.txt");
    fs3.open ("WFSin.txt");
    fs4.open ("WFCos.txt");

// %%%%%%%%%% writing data to files %%%%%%%%%%
// Phases + Probability distributions for all accepted collections of positions

    double Prob;
    double WF_RE, WF_IM, Phase;
    vector<double> XX;
    double PlotPoints = 30.; // number of plot points (in fact PlotPoints + 1 because of additional 0 in the loop)

    cout << "\n" << Chain.size() << "\n";
    cout << "Writing data to files" << "\n";

    for(unsigned int j = 0; j < Chain.size(); j++)
    {
        XX = Chain.at(j);

        for(double i = 0; i <= PlotPoints; i++ )
        {
        XX[0] = (i/PlotPoints)*L; // "positions" of the "last" particle
        WF_RE = _Permutator.WaveFunction(N, c, XX, K).real();
        WF_IM = _Permutator.WaveFunction(N, c, XX, K).imag();
        Prob = pow(WF_RE,2) + pow(WF_IM,2);

        if( WF_RE > 0){ // phase calculation
        Phase = atan( WF_IM/WF_RE );
        }
        else{
        Phase = atan( WF_IM/WF_RE ) + M_PI;
        }

        fs1 << Prob << " ";
        fs2 << Phase << " ";
        }

    fs1 << "\n";
    fs2 << "\n";
    cout << j << " '";
    }
    fs1.close();
    fs2.close();


*/





/*



//***********************************************************************************************************
//***********************************************************************************************************

// ########################## txt files preparation ##########################

    std::ofstream create1 ("WFS.txt");
    std::ofstream create2 ("WFPhase.txt");
    std::ofstream create3 ("Minima.txt");
    std::ofstream create4 ("WFSin.txt");
    std::ofstream create5 ("WFCos.txt");
    std::fstream fs1;
    std::fstream fs2;
    std::fstream fs3;
    std::fstream fs4;
    std::fstream fs5;
    fs1.open ("WFS.txt");
    fs2.open ("WFPhase.txt");
    fs3.open ("Minima.txt");
    fs4.open ("WFSin.txt");
    fs5.open ("WFCos.txt");

// ########################## writing data to files ##########################
// Phases, Probability distributions and Minima for all accepted collections of positions

    cout << "\n" << "\n" << "Phases/Probabilites/Minima calculation" << "\n"  ;

    int NumOfBisections = 8;
    int InitialDivision = 70;
    //int ChainElement = 20;
    double JumpRestriction = 0.7;

    cout << "\n" << Chain.size() << "\n";


    for(unsigned int i = 0; i < Chain.size(); i++ )
    {
    _Permutator.MinimaFinder(Minima,PhaseMatrix,WFMatrix,SinMatrix,CosMatrix,NumOfBisections,InitialDivision,N,i,JumpRestriction,c,L,Chain,K);

        if (i % 10 == 0)
        {
            cout << i << ", ";
        }

        // ======================= PROBABILITIES ======================
        for(int j = 0; j < InitialDivision; j++ )
        {
            fs1 << ( WFMatrix.back() ).at(j) << " ";
        }
        fs1 << "\n";
        // ============================================================
        // ========================== PHASES ==========================
        for(int j = 0; j < InitialDivision; j++ )
        {
            fs2 << ( PhaseMatrix.back() ).at(j) << " ";
        }
        fs2 << "\n";
        // ============================================================
        // ========================== MINIMA ==========================
        fs3 << Minima.back() << " ";
        // ============================================================
        // ========================== SIN =============================
        for(int j = 0; j < InitialDivision; j++ )
        {
            fs4 << ( SinMatrix.back() ).at(j) << " ";
        }
        fs4 << "\n";
        // ============================================================
        // ========================== COS =============================
        for(int j = 0; j < InitialDivision; j++ )
        {
            fs5 << ( CosMatrix.back() ).at(j) << " ";
        }
        fs5 << "\n";
        // ============================================================
    }

    fs1.close();
    fs2.close();
    fs3.close();
    fs4.close();
    fs5.close();


// ========== BEEP ===========
    cout << '\a';
// ===========================



*/




/* ################################################################################################################
***************************** SCALAR PRODUCTS CALCULATION and save to files ***************************************
################################################################################################################ */



/* ------------------------------------------------------------------------------------------------------------------------------
============== File with quasimomenta collections reading (K_Matrix) ============================================================
------------------------------------------------------------------------------------------------------------------------------ */


    std::ifstream k_matrix;
    k_matrix.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\K_3.txt");

    int NumOfQuasimomenta = 3;
    int NumOfQuasimomentaSets = 969;

    vector<vector<double>> K_Matrix;
    vector<double> TmpVec;
    double tmp;

    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        for(int i = 0; i < NumOfQuasimomenta; i++)
        {
             k_matrix >> tmp;
             TmpVec.push_back(tmp);
        }
        K_Matrix.push_back( TmpVec );
        TmpVec.clear();
    }

    k_matrix.close();

/*
    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        for(int i = 0; i < NumOfQuasimomenta; i++)
        { cout << K_Matrix[j][i] << ", " }
        cout << "\n";
    }
*/



/* ------------------------------------------------------------------------------------------------------------------------------
============== Scalar products calculation  =====================================================================================
------------------------------------------------------------------------------------------------------------------------------ */

    vector<vector<double>> AllX;
    vector<double> AllWs, AllNorms;

    int PolynomialDegree = 5;
    int IntervalDivision = 4;
    vector<double> K1 = K_Matrix.at(0);
    vector<double> K2 = K;
    vector<double> Measured_X = {.1,.2,.3};

//======= Preparation of collections which are indispensable to perform summations ===============================

    _Permutator.IntegrationCollections(K1,K2,AllX,AllWs,AllNorms,L,PolynomialDegree,IntervalDivision);

//----------------------------------------------------------------------------------------------------------------

    std::ofstream scalar_re ("Scalar_re_test.txt");
    std::ofstream scalar_im ("Scalar_im_test.txt");
    std::fstream sc_re;
    std::fstream sc_im;
    sc_re.open ("Scalar_re_test.txt");
    sc_im.open ("Scalar_im_test.txt");

    vector<double>  ScalarProducts_re, ScalarProducts_im, Sc_re_tmp, Sc_im_tmp;
    vector<int> Parallel_Order;
    Sc_re_tmp.reserve( NumOfQuasimomentaSets );
    Sc_im_tmp.reserve( NumOfQuasimomentaSets );
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Parallel for loop
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #pragma omp parallel for
    for(int i = 0; i < NumOfQuasimomentaSets; i++)
    {
        K1 = K_Matrix.at(i);
        double norm1 = Norms.at(i);
        double norm2 = norm;

        _Permutator.ScalarProduct(K1, K2, AllX, AllWs, AllNorms, c, Measured_X, ScalarProducts_re, ScalarProducts_im, norm1 , norm2);

        Parallel_Order.push_back(i);
        cout << " \n \n \n " << ScalarProducts_re.back() << " \n \n \n ";
        cout << "\n \n" << i << "\n \n";
    }
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Sort and writing to file:)
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        int position = Parallel_Order.at(j);

        Sc_re_tmp[ position ] = ScalarProducts_re.at( position );
        Sc_im_tmp[ position ] = ScalarProducts_im.at( position );
    }
    ScalarProducts_re = Sc_re_tmp;
    ScalarProducts_im = Sc_im_tmp;

    for(int s = 0; s < NumOfQuasimomentaSets; s++)
    {
        sc_re << ScalarProducts_re.at(s) << " ";
        sc_im << ScalarProducts_im.at(s) << " ";
    }


    //sc_re << ScalarProducts_re.back() << " ";
    //sc_im << ScalarProducts_im.back() << " ";
    //cout << " \n \n \n " << ScalarProducts_re.back() << " \n \n \n ";
    //cout << "\n \n" << i << "\n \n";

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    sc_re.close();
    sc_im.close();









/* ################################################################################################################
**************************************** TIME EVOLUTION ***********************************************************
################################################################################################################ */

//=================================================================================================================
//FIRSTLY WE NEED TO CHOOSE RELEVANT COEFFICIENTS AND ROWS OF K_MATRIX
//=================================================================================================================

//-----------------------------------------------------------------------------------------------------------------
//FILES READING
//-----------------------------------------------------------------------------------------------------------------
/*
    std::ifstream k_matrix;
    std::ifstream scprod_re;
    std::ifstream scprod_im;
    k_matrix.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\K_3.txt");
    scprod_re.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\Scalar_re_4.txt");
    scprod_im.open("C:\\Users\\Andrzej\\Desktop\\build-Permutator-Desktop_Qt_5_5_0_MinGW_32bit-Debug\\Scalar_im_4.txt");

    int NumOfQuasimomenta = 3;
    int NumOfQuasimomentaSets = 969;

    vector<vector<double>> K_Matrix;
    vector<double> TmpVec, ScalarProducts_re, ScalarProducts_im;
    double tmp_k, tmp_sc_re, tmp_sc_im;

    for(int j = 0; j < NumOfQuasimomentaSets; j++)
    {
        for(int i = 0; i < NumOfQuasimomenta; i++)
        {
             k_matrix >> tmp_k;
             TmpVec.push_back(tmp_k);
        }
        K_Matrix.push_back( TmpVec );
        TmpVec.clear();

        scprod_re >> tmp_sc_re;
        scprod_im >> tmp_sc_im;

        ScalarProducts_re.push_back(tmp_sc_re);
        ScalarProducts_im.push_back(tmp_sc_im);
    }

    k_matrix.close();
    scprod_re.close();
    scprod_im.close();




//-----------------------------------------------------------------------------------------------------------------
// WE CHOOSE RELEVANT ELEMENTS
//-----------------------------------------------------------------------------------------------------------------

    double CutOff_Coefficient =  0.0002;

    _Permutator.RelevantCoefficients(ScalarProducts_re, ScalarProducts_im, K_Matrix, CutOff_Coefficient, Norms);

    cout << "Number of relevant elements: " << ScalarProducts_im.size() << "\n";


//-----------------------------------------------------------------------------------------------------------------
// METROPOLIS FOR DIFFERENT TIME MOMENTS
//-----------------------------------------------------------------------------------------------------------------


    int Steps_t = 50000;

    int N_Remaining = (K_Matrix.at(0)).size();

    vector<vector<vector<double>>> Chain_T_dependence;


    std::ofstream chainfile ("X_T.txt");
    std::fstream chain_t;
    chain_t.open ("X_T.txt");

    vector<vector<double>> Chain_t;

    for(int j = 0; j < 1; j++)
    {
        double t = 0.7 * j;

        cout << "\n\n" << "time value: " << t << "\n\n";

        //vector<vector<double>> Chain_t;
        vector<double> WaveFValues_t;
        vector<double> X0_t = {.012, .95, .33};
        vector<double> X1_t = X0_t;
        _Permutator.AddCollection(Chain_t, X0_t);

        for(int i = 0; i < Steps_t; i++)
        {

        _Permutator.Metropolis_TimeDependent(N_Remaining, c, t, Chain_t, WaveFValues_t, X0_t, X1_t, K_Matrix, ScalarProducts_re, ScalarProducts_im, delta, L, Norms);

            if (i % 100 == 0)
            {
                cout << i << ", ";
            }
        }


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

    int NumOfPoints = 100;
    double WF6, WF3;
    double Wf6Value = 0.;
    double Wf3Value = 0.;


    std::ofstream wf6 ("wf6.txt");
    std::ofstream wf3 ("wf3.txt");

    std::fstream w6;
    std::fstream w3;

    w6.open ("wf6.txt");
    w3.open ("wf3.txt");



    for(int j = 0; j < NumOfPoints ; j++)
    {
        vector<double> X_01 = {.1, .2, .3, .12, .75 };
        vector<double> X_02 = {.12, .75 };

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


















/* ================================================================================================================
###################################################################################################################
----------------------------- REPULSIVE CASE ----------------------------------------------------------------------
###################################################################################################################
================================================================================================================ */

/*



    vector<double> K_re = {0.000150252, 0.0000977129, 0.0000552789, 0.000017858, -0.0000179714, -0.0000553577, -0.0000976897, -0.000150083};
    vector<double> K_im = {-4.37684, -2.97843, -1.74707, -0.576849, 0.57685, 1.74707, 2.97843, 4.37684};
    vector<double> R_X = {.1, .2, .3, .4, .5, .6, .7, .8};

    std::complex<double> R_Wave;

    R_Wave = _Permutator.R_WaveFunction(N , c, R_X, K_re, K_im);

    cout<< R_Wave.real() << ", " << R_Wave.imag();


    vector<vector<double>> AllX;
    vector<double> AllWs, AllNorms;

    int PolynomialDegree = 10;
    int IntervalDivision = 5;
    vector<double> K1_re = {.5, .3, 9.};
    vector<double> K1_im = {0, .45, 19.};
    vector<double> K2_re = {1. ,2. ,3. ,0. ,3. ,12. };
    vector<double> K2_im = {1. , 1.,.5 ,2. ,7. , 0.};
    vector<double> Measured_X = {.2,.4,.3};

    _Permutator.R_IntegrationCollections(K1_re, K2_re ,AllX ,AllWs ,AllNorms ,L, PolynomialDegree,  IntervalDivision);

    _Permutator.R_ScalarProduct(K1_re, K1_im, K2_re, K2_im, AllX, AllWs, AllNorms, c, Measured_X);


*/





double t1=clock();
cout<< "\n\n\n"<< (t1-t0)/CLOCKS_PER_SEC << "\n\n";

}



