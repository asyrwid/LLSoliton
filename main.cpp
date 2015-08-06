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
    srand( time( NULL ) );

/* ************** INITIAL VALUES AND OBJECTS ************** */
//***********************************************************************************************************

    int    Steps = 100000; // number of metropolis steps
    int    N = 8; // number of particles
    double c = 16.0; // coupling constant
    double L = 1.0; // system length
    double delta = 0.2; // maximal "jump" between single particle position from one realization to another

    // initial collection of particle positions
    vector<double> X0 = {0.5,.1,.15,.2,.25,.3,.35,.4};
    vector<double> X1 = X0;

    // quasimomentas
    vector<double> K = {-16.6144, -12.5642, -8.62718, -4.72427, 11.0075, 14.9104, 18.8474, 22.8976};
    // chain of position
    vector<vector<double> > Chain;

    //chain of wave function values
    vector<double> WaveFValues;

    //vector of wave function minima and Phase/WF PlotPoints
    vector<double> Minima;
    vector<vector<double> > PhaseMatrix;
    vector<vector<double> > WFMatrix;

//***********************************************************************************************************

permutator _Permutator;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //cout << _Permutator.getPermSize(); // initially _N=0;
    _Permutator.setPermSize(N); // de facto particle number
    //cout << _Permutator.getPermSize(); // N is established now!
    _Permutator.Permute(N,N); // all permutations preparation
    _Permutator.AddCollection(Chain, X0);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// ************** METROPOLIS **************
//***********************************************************************************************************

    for(int i = 0; i < Steps; i++)
    {
    _Permutator.Metropolis(N,c,Chain,WaveFValues,X0,X1,K,delta,N,L);
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







/*



// ****** LAST PROBABILITY AND PHASE OF WAVE FUNCTION ******
//***********************************************************************************************************

// %%%%%%%%%% txt files preparation %%%%%%%%%%
    std::ofstream create1 ("WFS.txt");
    std::ofstream create2 ("WFPhase.txt");
    std::fstream fs1;
    std::fstream fs2;
    fs1.open ("WFS.txt");
    fs2.open ("WFPhase.txt");

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







//***********************************************************************************************************
//***********************************************************************************************************

// ########################## txt files preparation ##########################

    std::ofstream create1 ("WFS.txt");
    std::ofstream create2 ("WFPhase.txt");
    std::ofstream create3 ("Minima.txt");
    std::fstream fs1;
    std::fstream fs2;
    std::fstream fs3;
    fs1.open ("WFS.txt");
    fs2.open ("WFPhase.txt");
    fs3.open ("Minima.txt");

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
    _Permutator.MinimaFinder(Minima,PhaseMatrix,WFMatrix,NumOfBisections,InitialDivision,N,  i  ,JumpRestriction,c,L,Chain,K);

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

    }


    fs1.close();
    fs2.close();
    fs3.close();



/*


// ########################## txt files preparation ##########################

    std::ofstream create1 ("WFS.txt");
    std::ofstream create2 ("WFPhase.txt");
    std::ofstream create3 ("Minima.txt");
    std::fstream fs1;
    std::fstream fs2;
    std::fstream fs3;
    fs1.open ("WFS.txt");
    fs2.open ("WFPhase.txt");
    fs3.open ("Minima.txt");

// ########################## writing data to files ##########################
// Phases, Probability distributions and Minima for all accepted collections of positions


// ======================+ PROBABILITIES ======================
    for(unsigned i = 0; i < WFMatrix.size(); i++)
    {
        for(int j = 0; j < InitialDivision; j++ )
        {
            fs1 << (WFMatrix.at(i)).at(j) << " ";
        }
        fs1 << "\n";
    }
    fs1.close();


// ========================== PHASES ==========================
    for(unsigned i = 0; i < PhaseMatrix.size(); i++)
    {
        for(int j = 0; j < InitialDivision; j++ )
        {
            fs2 << (PhaseMatrix.at(i)).at(j) << " ";
        }
        fs2 << "\n";
    }
    fs2.close();



// ========================== MINIMA ==========================
    for(unsigned int i = 0; i < Minima.size(); i++)
    {
            fs3 << Minima.at(i) << " ";
    }
    fs3 << "\n";
    fs3.close();


// ========== BEEP ==========
    cout << '\a';
// ==========================





*/


}



