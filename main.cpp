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

/* ************** INITIAL VALUES ************** */
//***********************************************************************************************************

    int    Steps = 200; // number of metropolis steps
    int N = 8; // number of particles
    double c = 0.08; // coupling constant
    double L = 1; // system length
    double delta = 0.3; // maximal "jump" between single particle position from one realization to another

    // initial collection of particle positions
    vector<double> X0 = {.1,.12,.13,.122,.132,.111,.113,.118};
    vector<double> X1 = X0;

    // quasimomentas
    vector<double> K = {-0.752974,-0.30724,0.107856,0.552964,5.73022,6.17533,6.59043,7.03616};

    // chain of position
    vector<vector<double> > Chain;

    //chain of wave function values
    vector<double> WaveFValues;


    vector<double> Jumps;

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
    cout << i << ", ";
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
    double Prob;
    double WF_RE, WF_IM, Phase;
    for(double i = 0; i <= 100; i++ )
    {
    X0[0] = (i/100.)*L; // "positions" of the "last" particle
    WF_RE = _Permutator.WaveFunction(N, c, X0, K).real();
    WF_IM = _Permutator.WaveFunction(N, c, X0, K).imag();
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
    fs1.close();
    fs2.close();

//***********************************************************************************************************

    cout << "\n"  << "\n"  << "\n"  ;
    _Permutator.PhaseJump(Jumps, 10, N, 10,1.0 ,c,L,Chain,K);
    cout << "\n"  << "\n"  << "\n"  ;
    for(unsigned i = 0; i < Jumps.size(); i++)
    {
    cout << Jumps.at(i) << "\n"  ;
    }








}



