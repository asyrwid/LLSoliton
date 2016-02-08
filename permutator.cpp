#include "permutator.h"
#include <iterator>
#include <math.h>
#include <random>
#include <omp.h>
#include <iostream>
#include <fstream>

//permutator::permutator()
//{
//}


# define Pi  3.14159265358979323846




void permutator::SaveVecToFile(vector<double> from,
                               std::ofstream &to)
{
    int veclength = from.size();
    for(int i = 0; i < veclength; i++)
    {
       to << from[i] << " ";
    }
}



void permutator::SaveVecFromFile(int veclength,
                                 std::ifstream &from,
                                 vector<double> &to)
{
    for(int i = 0; i < veclength; i++)
    {
       double tmp;
       from >> tmp;
       to.push_back(tmp);
    }
}



void permutator::SaveMatrixFromFile(int rows,
                                    int columns,
                                    std::ifstream &from,
                                    vector<vector<double>> &to)
{
    for(int i = 0; i < rows; i++)
    {
       vector<double> tmpvec;

       for( int j = 0; j < columns; j++)
       {
           double tmp;
           from >> tmp;
           tmpvec.push_back(tmp);
       }
       to.push_back(tmpvec);
    }
}



void permutator::SaveMatrixToFile(vector<vector<double>> from,
                             std::ofstream &to)
{
    int rows = from.size();
    int columns = (from.at(0)).size();

    for(int i = 0; i < rows; i++)
    {
       for(int j = 0; j < columns; j++)
       {
          to << (from.at(i)).at(j) << " ";
       }
       to << "\n";
    }
}













void permutator::setPermSize(int N)   //particle number setting
{                                     //and initial set 1,...,N
    _N = N;                           //preparation
    _PermPool.clear();
    for (int i=0; i<N; i++)
    {
        _PermPool.push_back(i);
    }
}




int permutator::getPermSize() //expectorate (wykrztuś,wypluj)
{                             //the particle number
    return _N;
}




void permutator::PermSwap(int a, int b) //"swap" additional operator
{
    if (a == b)
    {
     return;
    }

//    _PermPool.swap(a,b);
    std::swap(_PermPool[a],_PermPool[b]);
    _Permity = !_Permity; //de facto Parity change after

    //swap operation :)
}



void permutator::appendPerm() //operation appending
{                                //another permutation + sign
    //for(int i=0;i<size;i++)
        Perms.push_back(_PermPool);
    //cout << _PermPool;
    if (_Permity == true){
        _PermSigns.push_back(1);
    }
    else{
        _PermSigns.push_back(-1);
    }
    //cout << _PermSigns.last();
    //cout << Perms.last() << ' ' << _PermSigns.last();
}



void permutator::Permute(int k, int size)//all permutations
{                        //preparation + signs (push_backPerm used)
    if(k == 0)
    {
        appendPerm();
    }
    else{
       for(int i = k-1; i >= 0; i--)
            {
            PermSwap(i,k-1);
            Permute(k-1,size);
            PermSwap(i,k-1);
            }
    }
}



int permutator::GetPermEl(int i, int j) // acces to (i,j) element of
{                                        // all permutations array
    return (Perms.at(i)).at(j);
}




int permutator::SignFunction(double x1, double x2) // sign function
{                                                   // sign(x1-x2)
    int t;
    if(x1 == x2)
    {
        t = 0;
    }
    else
    {
    if(x1 < x2){
        t = -1;}
    else{
        t = 1;}
    }
    return t;
}




//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//*********************************** WAVE FUNCTION PREPARATION **********************************************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

typedef std::complex<double> Complex;



std::complex<double> permutator::sProduct(int n,
                                          int PermNumber,
                                          double c,
                                          vector<double> X,
                                          vector<double> K)
{
// ====================================================================================================================================
// 1) Product (j>k which appears in the wave function) definition. PermNumber is the number of permutation in permutations array (row)
// ====================================================================================================================================

    Complex  prod(1.,0.);
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < j; k++)
        {
        prod = prod* Complex( K[GetPermEl(PermNumber,j)] - K[GetPermEl(PermNumber,k)], -c*SignFunction( X[j], X[k]) ) ;
        }

    }
    // cout << prod.real() << ',' << prod.imag();


// ====================================================================================================================================
// 2) Exponent with a sum over X[]K[]
// ====================================================================================================================================

    double sumexp = 0.;
    for(int j = 0; j < n; j++)
    {
        sumexp = sumexp + K[GetPermEl(PermNumber,j)] * X[j];
    }
    Complex ImSumExp(0., sumexp);
    Complex ImExp;
    ImExp = exp(ImSumExp);
    // cout << ImExp.real() << ',' << ImExp.imag();

// ====================================================================================================================================
// 3) Sign of the permutation given by PermNumber
// ====================================================================================================================================

    int SPV = _PermSigns.at(PermNumber); // + or - 1 (Parity)
    Complex SignValue(SPV, 0.);

// ====================================================================================================================================
// 4) Finally, we can calculate one element of "big" sum apearing in the wave function
// ====================================================================================================================================

    Complex ProductElement;
    ProductElement = SignValue * ImExp * prod;

/*
    cout <<'('<<ProductElement.real()<<','<<ProductElement.imag()<<')' ;
    cout <<'('<<SignValue.real()<<','<<SignValue.imag()<<')' ;
    cout << prod.real() << ImExp.real();
    cout << prod.imag() << ImExp.imag();
*/

   // cout << ProductElement.real() << ',' << ProductElement.imag();
    return ProductElement;
}






std::complex<double> permutator::sProduct_unnormalized(int n,
                                          int PermNumber,
                                          double c,
                                          vector<double> X,
                                          vector<double> K)

//sProduct divided by c (if c is big - we should do it cause there is problem with normalization)
{
// ====================================================================================================================================
// 1) Product (j>k which appears in the wave function) definition. PermNumber is the number of permutation in permutations array (row)
// ====================================================================================================================================

    Complex  prod(1.,0.);
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < j; k++)
        {
        prod = prod* Complex( (K[GetPermEl(PermNumber,j)] - K[GetPermEl(PermNumber,k)])/c, -SignFunction( X[j], X[k]) ) ;
        }

    }
    // cout << prod.real() << ',' << prod.imag();


// ====================================================================================================================================
// 2) Exponent with a sum over X[]K[]
// ====================================================================================================================================

    double sumexp = 0.;
    for(int j = 0; j < n; j++)
    {
        sumexp = sumexp + K[GetPermEl(PermNumber,j)] * X[j];
    }
    Complex ImSumExp(0., sumexp);
    Complex ImExp;
    ImExp = exp(ImSumExp);
    // cout << ImExp.real() << ',' << ImExp.imag();

// ====================================================================================================================================
// 3) Sign of the permutation given by PermNumber
// ====================================================================================================================================

    int SPV = _PermSigns.at(PermNumber); // + or - 1 (Parity)
    Complex SignValue(SPV, 0.);

// ====================================================================================================================================
// 4) Finally, we can calculate one element of "big" sum apearing in the wave function
// ====================================================================================================================================

    Complex ProductElement;
    ProductElement = SignValue * ImExp * prod;

/*
    cout <<'('<<ProductElement.real()<<','<<ProductElement.imag()<<')' ;
    cout <<'('<<SignValue.real()<<','<<SignValue.imag()<<')' ;
    cout << prod.real() << ImExp.real();
    cout << prod.imag() << ImExp.imag();
*/

   // cout << ProductElement.real() << ',' << ProductElement.imag();
    return ProductElement;
}







std::complex<double> permutator::sProduct_strong(int n,
                                                 int PermNumber,
                                                 double c,
                                                 vector<double> X,
                                                 vector<double> K)
// strong sProduct is a function which I use in the case of very strong coupling, where I assume that c^2 > > (k_j1-k_s1)(k_j2-k_s2)
// therefore I will consider only two terms from the full product - proporional to c^N(N-1)/2 and to c^[N(N-1)/2-1]
// It should work very well cause in the c->oo limit k-> 2pi/L * integer
// We divide our wave function by c^N(N-1)/2 because... we can :) It is just a constant factor which can be troublesome (range of double...)
// prof. Gajda hint

{

// ====================================================================================================================================
// 1) Product (j>k which appears in the wave function) definition. PermNumber is the number of permutation in permutations array (row)
// ====================================================================================================================================

    Complex  prodcc(1.,0.);
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < j; k++)
        {
        prodcc = prodcc* Complex( 0. , - SignFunction( X[j], X[k]) ) ;
        }
    }


    Complex  prodck(0.,0.);
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < j; k++)
        {
        prodck = prodck + prodcc/ Complex(0., -c*SignFunction(X[j], X[k]) ) * Complex( K[GetPermEl(PermNumber,j)] - K[GetPermEl(PermNumber,k)] , 0. ) ;
        }
    }


    Complex  prod;

    prod = prodcc + prodck;

// ====================================================================================================================================
// 2) Exponent with a sum over X[]K[]
// ====================================================================================================================================

    double sumexp = 0.;
    for(int j = 0; j < n; j++)
    {
        sumexp = sumexp + K[GetPermEl(PermNumber,j)] * X[j];
    }
    Complex ImSumExp(0., sumexp);
    Complex ImExp;
    ImExp = exp(ImSumExp);
    // cout << ImExp.real() << ',' << ImExp.imag();

// ====================================================================================================================================
// 3) Sign of the permutation given by PermNumber
// ====================================================================================================================================

    int SPV = _PermSigns.at(PermNumber); // + or - 1 (Parity)
    Complex SignValue(SPV, 0.);

// ====================================================================================================================================
// 4) Finally, we can calculate one element of "big" sum apearing in the wave function
// ====================================================================================================================================

    Complex ProductElement;
    ProductElement = SignValue * ImExp * prod;

/*
    cout <<'('<<ProductElement.real()<<','<<ProductElement.imag()<<')' ;
    cout <<'('<<SignValue.real()<<','<<SignValue.imag()<<')' ;
    cout << prod.real() << ImExp.real();
    cout << prod.imag() << ImExp.imag();
*/

   // cout << ProductElement.real() << ',' << ProductElement.imag();
    return ProductElement;
}










std::complex<double> permutator::WaveFunction(int n,
                                              double c,
                                              vector<double> X,
                                              vector<double> K,
                                              double norm)
{
    // N! calculation - we need to know how many permutations we will have
    //*************************************************************************
            int Nfact = 1; // Nfact = N! (Nfactorial)
            for (int i = 2; i <= n; i++) // Nfactorial value (N!) calculation
            { Nfact = Nfact * i; }
            //cout << Nfact;
    // Actually, it is the number of rows of "perms"
    // ************************************************************************

    double re = 0., im = 0.;
    int i;

    //omp_set_dynamic(1);
    //omp_set_num_threads(2);
    //#pragma omp parallel private(i) reduction(+:re) reduction(+:im)
    //{
        Complex wavefunction(0.,0.);
        //#pragma omp for
        for(i = 0; i < Nfact; i++)
        {
            wavefunction = wavefunction + sProduct(n, i, c, X, K);
            //cout<< i << "\n" ;
        }
        re=wavefunction.real();
        im=wavefunction.imag();
    //}
    Complex Wavefunction(re,im);

    // Normalization factor

    double norm2 = Nfact;
    for (int j = 0; j < n ; j++)
    {
        for (int k = 0; k < j; k++)
        {
        norm2 = norm2 * ( pow ((K[j]-K[k]),2.) + pow (c,2.) );
        }
    }


    Complex Norm( pow(norm2, - 0.5) * pow(norm, -0.5) ,0. ) ;
    Wavefunction = Norm * Wavefunction;


    //cout << fabs(wavefunction) <<  "\n";

    //cout << wavefunction.real() << ',' << wavefunction.imag();

    //wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return Wavefunction;
}




std::complex<double> permutator::WaveFunction_unnormalized(int n,
                                              double c,
                                              vector<double> X,
                                              vector<double> K,
                                              double norm)
{
    // N! calculation - we need to know how many permutations we will have
    //*************************************************************************
            int Nfact = 1; // Nfact = N! (Nfactorial)
            for (int i = 2; i <= n; i++) // Nfactorial value (N!) calculation
            { Nfact = Nfact * i; }
            //cout << Nfact;
    // Actually, it is the number of rows of "perms"
    // ************************************************************************

    double re = 0., im = 0.;
    int i;

    //omp_set_dynamic(1);
    //omp_set_num_threads(2);
    //#pragma omp parallel private(i) reduction(+:re) reduction(+:im)
    //{
        Complex wavefunction(0.,0.);
        //#pragma omp for
        for(i = 0; i < Nfact; i++)
        {
            wavefunction = wavefunction + sProduct_unnormalized(n, i, c, X, K);
            //cout<< i << "\n" ;
        }
        re=wavefunction.real();
        im=wavefunction.imag();
    //}
    Complex Wavefunction(re,im);

    // Normalization factor
/*
    double norm2 = Nfact;
    for (int j = 0; j < n ; j++)
    {
        for (int k = 0; k < j; k++)
        {
        norm2 = norm2 * ( pow ((K[j]-K[k]),2.) + pow (c,2.) );
        }
    }


    Complex Norm( pow(norm2, - 0.5) * pow(norm, -0.5) ,0. ) ;
    Wavefunction = Norm * Wavefunction;
*/

    Complex Norm( /*pow(norm2, - 0.5)*/  pow(norm, -0.5) ,0. ) ;
    Wavefunction = Norm * Wavefunction;

    //cout << fabs(wavefunction) <<  "\n";

    //cout << wavefunction.real() << ',' << wavefunction.imag();

    //wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return Wavefunction;
}



std::complex<double> permutator::WaveFunction_strong(int n,
                                              double c,
                                              vector<double> X,
                                              vector<double> K,
                                              double norm)
{
    // N! calculation - we need to know how many permutations we will have
    //*************************************************************************

            int Nfact = 1; // Nfact = N! (Nfactorial)
            for (int i = 2; i <= n; i++) // Nfactorial value (N!) calculation
            { Nfact = Nfact * i; }

            //cout << Nfact;
    // Actually, it is the number of rows of "perms"
    // ************************************************************************

    double re = 0., im = 0.;
    int i;

    //omp_set_dynamic(1);
    //omp_set_num_threads(2);
    //#pragma omp parallel private(i) reduction(+:re) reduction(+:im)
    //{
        Complex wavefunction(0.,0.);
        //#pragma omp for
        for(i = 0; i < Nfact; i++)
        {
            wavefunction = wavefunction + sProduct_strong(n, i, c, X, K);
            //cout<< i << "\n" ;
        }
        re=wavefunction.real();
        im=wavefunction.imag();
    //}
    Complex Wavefunction(re,im);

    // Normalization factor
/*
    double norm2 = Nfact;
    for (int j = 0; j < n ; j++)
    {
        for (int k = 0; k < j; k++)
        {
        norm2 = norm2 * ( pow ((K[j]-K[k]),2.) + pow (c,2.) );
        }
    }
*/

    Complex Norm( /*pow(norm2, - 0.5)*/  pow(norm, -0.5) ,0. ) ;
    Wavefunction = Norm * Wavefunction;


    //cout << fabs(wavefunction) <<  "\n";

    //cout << wavefunction.real() << ',' << wavefunction.imag();

    //wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return Wavefunction;
}








//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//*********************************** COLLECTIONS OF PARTICLE POSITIONS GENERATION ***************************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

double permutator::GetRandom(double min, double max)
{// Returns a random double between min and max
    return ((double) rand()*(max-min)/(double)RAND_MAX + min);
}





void permutator::Positions(vector<double> X0,
                           vector<double> &X1,
                           double delta,
                           int n,
                           double L)
{// one step in Metropolis x1(n)->x1'(n) = x0(n)+ delta*Q, Q in [-L,L]
    double posX;
        posX = delta * GetRandom(-L, L);
        if( X0.at(n) + posX < 0 ){
        X1[n] = L + ( X0.at(n) + posX);
        }
        else
        {
            if( X0.at(n) + posX > L )
            {
            X1[n] =  X0.at(n) + posX - L;
            }
            else
            {
            X1[n] =  X0.at(n) + posX;
            }
        }
    //cout << X1;
}





void permutator::AddCollection(vector<vector<double> > &Chain,
                               vector<double>  Set)
{
    Chain.push_back( Set );
    //cout << Chain;
}





void permutator::Metropolis(int steps,
                            int N,
                            double c,
                            vector<vector<double> > &Chain,
                            vector<double> &WaveFValues,
                            vector<double> &X0,
                            vector<double> &X1,
                            vector<double> K,
                            double delta,
                            double L,
                            double norm)
{
    double ratio;
    double dif;
    double random0to1;
    double wf0;
    double wf1;

    wf0 =  pow (fabs(WaveFunction(N,c, X0, K, norm) ), 2.);

for(int u = 0; u < steps; u++)
{


    for(int i = 0; i < N; i++)
    {
        Positions(X0, X1, delta, i, L);
    }
        wf1 =  pow (fabs(WaveFunction(N,c, X1, K, norm) ), 2.);

        ratio = wf1/wf0;

        if( ratio >= 1. )
        {
            WaveFValues.push_back(wf1);
            AddCollection(Chain,X1);
        }
        else
        {
            random0to1 = GetRandom(0.,1.);
            dif = ratio - random0to1;
                if( dif > 0. )
                {
                    WaveFValues.push_back(wf1);
                    AddCollection(Chain,X1);
                }
                else
                {
                    X1 = X0;
                    wf1 = wf0;
                }
        }
     X0 = X1;
     wf0 = wf1;

}

}





void permutator::Metropolis_unnormalized(int OutSteps,
                                   int InsideSteps,
                            int N,
                            double c,
                            vector<vector<double> > &Chain,
                            vector<double> &WaveFValues,
                            vector<double> &X1,
                            vector<double> K,
                            double delta,
                            double L,
                            double norm)
{
    double ratio;
    double dif;
    double random0to1;
    double wf0;
    double wf1;

for(int f = 0; f < OutSteps; f++)
{
    vector<double> X0;

    for(int z = 0; z < N; z++)
    {
        X0.push_back(GetRandom(0.,1.));
    }

    wf0 =  pow (fabs(WaveFunction_unnormalized(N,c, X0, K, norm) ), 2.);

for(int u = 0; u < InsideSteps; u++)
{


    for(int i = 0; i < N; i++)
    {
        Positions(X0, X1, delta, i, L);
    }
        wf1 =  pow (fabs(WaveFunction_unnormalized(N,c, X1, K, norm) ), 2.);

        ratio = wf1/wf0;

        if( ratio >= 1. )
        {
            WaveFValues.push_back(wf1);
            AddCollection(Chain,X1);
            X0 = X1;
            wf0 = wf1;
        }
        else
        {
            random0to1 = GetRandom(0.,1.);
            dif = ratio - random0to1;
                if( dif > 0. )
                {
                    WaveFValues.push_back(wf1);
                    AddCollection(Chain,X1);
                    X0 = X1;
                    wf0 = wf1;
                }
                else
                {
                    X1 = X0;
                    wf1 = wf0;
                }
        }

     if(u%100 == 0)
     {
        //cout<< "Out" << "\r"<< f << " ";
        cout<< "\r"<< u << " ";
     }
}

}


}






void permutator::Metropolis_strong(int OutSteps,
                                   int InsideSteps,
                            int N,
                            double c,
                            vector<vector<double> > &Chain,
                            vector<double> &WaveFValues,
                            vector<double> &X1,
                            vector<double> K,
                            double delta,
                            double L,
                            double norm)
{
    double ratio;
    double dif;
    double random0to1;
    double wf0;
    double wf1;

for(int f = 0; f < OutSteps; f++)
{
    vector<double> X0;

    for(int z = 0; z < N; z++)
    {
        X0.push_back(GetRandom(0.,1.));
    }

    wf0 =  pow (fabs(WaveFunction_strong(N,c, X0, K, norm) ), 2.);

for(int u = 0; u < InsideSteps; u++)
{


    for(int i = 0; i < N; i++)
    {
        Positions(X0, X1, delta, i, L);
    }
        wf1 =  pow (fabs(WaveFunction_strong(N,c, X1, K, norm) ), 2.);

        ratio = wf1/wf0;

        if( ratio >= 1. )
        {
            WaveFValues.push_back(wf1);
            AddCollection(Chain,X1);
            X0 = X1;
            wf0 = wf1;
        }
        else
        {
            random0to1 = GetRandom(0.,1.);
            dif = ratio - random0to1;
                if( dif > 0. )
                {
                    WaveFValues.push_back(wf1);
                    AddCollection(Chain,X1);
                    X0 = X1;
                    wf0 = wf1;
                }
                else
                {
                    X1 = X0;
                    wf1 = wf0;
                }
        }

     if(u%100 == 0)
     {
        //cout<< "Out" << "\r"<< f << " ";
        cout<< "\r"<< u << " ";
     }
}

}


}
























int permutator::MaximumPosition(vector<double> WaveFValues)
{//Finding the set of positions (obtained by Metropolis algorithm) which maximizes the probability
 //in order to prepare plots of probability denstity and phase distribution of last particle
    int i = 0;
    int length = WaveFValues.size();

    for(int j = 1; j < length ; j++)
    {
        if( WaveFValues.at(j) <= WaveFValues.at(i) )
        {}
        else
        {
            i = j;
        }
    }
    return i;
}





//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//***************************** FINDING THE MINIMUM OF "THE LAST" PROBABILITY DENSITY ************************************************
//**************************** WE NEED TO FIND A JUMP(S) OF THE PHASE OF WAVE FUNCTION ***********************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

vector<double> permutator::PhaseCalculation(std::complex<double> WF)
{// gives a vector (fabs(Wavefunction), Phase), value of Wave Function is needed to find appropriate jump
 // we can use collected values to plot Prob density and phase distribution
 // additionally we calculate sin and cos of phase

    double  Phase, Sin, Cos;
    vector<double> Out;

    Phase = arg ( WF );
    Sin = WF.imag() /( fabs( WF ) );
    Cos = WF.real() /( fabs( WF ) );

    Out.push_back( pow ( fabs (WF), 2. ) ); // Prob density value (position 0)
    Out.push_back( Phase ); // Phase value (position 1)
    Out.push_back( Sin ); // Sin value (position 2)
    Out.push_back( Cos ); // Cos value (position 3)

    return Out;
}





void permutator::MinimaFinder(vector<double> &Minima,
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
                              double norm)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The MinimaFinder function find the value of proper minimum of prob density based on the phase jump and value of prob
density in the point corresponding to the phase discontinuity. We know that the jump phase shouldn't be bigger than Pi :)

    Finding the values of phase and prob density in the InitialDivision number of points we can find an interesting
bracket. After that we use the bisection method. Phase and Prob density plot points obtained during initial calculation
may be used to prepare plots - so we can kill two birds with one stone.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

// ======================================================================================================================
// 1) Calculation of the Phase and Prob density values in initial points
// ======================================================================================================================

    vector<double>  PosX = Chain.at(ChainElement); // Choosing the collection of positions from the Chain (  .at(ChainElement)  )
    vector<double>  InitialPositions, PhaseValues, WFValues, SinValues, CosValues, temp;
    double  Phase, WF, SinValue, CosValue;

//Preparation of initial set of points (0,L*1/(InitialDiv-1),...,L)

    for(int i = 0; i < InitialDivision; i++)
    {
        InitialPositions.push_back( i*L/(InitialDivision - 1) );
    }

//Values of Phase and fabs(Wave Function)^2 in all the points from initial set

    for(int i = 0; i < InitialDivision; i++)
    {
        PosX.at(0) = InitialPositions.at(i);
        temp = PhaseCalculation(WaveFunction(n, c, PosX, K, norm) );

        WF = temp.at(0);
        Phase = temp.at(1);
        SinValue = temp.at(2);
        CosValue = temp.at(3);

        WFValues.push_back( WF );        // Set of Prob denisty values in initial points
        PhaseValues.push_back( Phase );  // Set of Phase values in initial points
        SinValues.push_back( SinValue ); // Set of Sin values in initial points
        CosValues.push_back( CosValue ); // Set of Cos values in initial points

       // cout << i << "   " << PhaseValues.at(i) << "   " << WFValues.at(i) <<"\n";
    }

//Now we want to append values of Phase and Prob (and Sin and Cos) density to matrices

    AddCollection(PhaseMatrix, PhaseValues);
    AddCollection(WFMatrix, WFValues);

    AddCollection(SinMatrix, SinValues);
    AddCollection(CosMatrix, CosValues);


// ======================================================================================================================
// 2) Proper value of the Phase Jump searching
// ======================================================================================================================

//Calculation of diferences between all 2 phase values: fabs( Phase(j+1)-Phase(j) ); and finding JUMPS!

    vector<double> JumpValue;
    vector<int> JumpPositions, JumpPositions2;
    double diff;

    for(int i = 0; i < InitialDivision-1; i++)
    {
        //cout << PhaseValues.at(i+1) << "   " << PhaseValues.at(i) <<  "   " <<   PhaseValues.at(i + 1) - PhaseValues.at( i )  <<"\n";
        diff =  fabs ( PhaseValues.at(i + 1) - PhaseValues.at( i ) );
        JumpValue.push_back( diff );

        if( diff >= JumpRestriction)
        {
            JumpPositions.push_back(i);
        }
    }

/*
    cout << "\n" << "-----1------" << "\n";

    for(unsigned i = 0; i < JumpPositions.size(); i ++)
    {
    cout << JumpPositions.at(i) << "\n";
    }
*/


//We need to find the smallest value of Prob density in the set of points around the JUMPS

    int s, p;
    for(unsigned int j = 0; j < JumpPositions.size(); j++)
    {
        s = JumpPositions.at(j);
        JumpPositions2.push_back(s);
        JumpPositions2.push_back(s+1);
    }

/*
    cout << "-----2------" << "\n";
    for(unsigned i = 0; i < JumpPositions2.size(); i++)
    {
        cout << JumpPositions2.at(i) << "\n";
    }
*/

    // p is the number of element corresponding to minimal value
    // of wave function inside the bracket corresponding to the jump of phase

    p = JumpPositions2.at(0);
    for(int j = 1; j < JumpPositions2.size(); j++)
    {
        s = JumpPositions2.at(j);
        if( WFValues.at(p) < WFValues.at(s) )
        {
        }
        else{
            p = s;
        }
    }
    //cout << "\n" << "============" << "\n";
    //cout << p << "\n";


// ======================================================================================================================
// 3) Last step - BISECTION METHOD
// ======================================================================================================================

// (A=p-1,B=p,C=p+1) st. f(A)>f(B),f(C)>f(B) - bisection method

    double A,B,C,D,E;
    int pa, pb = p, pc; // number of collection element corresponding to position A/B/C

    if ( p > 0 )
    {
        A = InitialPositions.at(p-1);
        pa = p - 1;
        if ( p < InitialDivision-1 )
            {
                C = InitialPositions.at(p+1);
                pc = p + 1;
            }
            else //Means p = InitialDivision-1
            {
                C = InitialPositions.at(0);
                pc = 0;
            }
    }else//Means p = 0
    {
        A = InitialPositions.at(InitialDivision-1);
        pa = InitialDivision - 1;
        C = InitialPositions.at(p+1);
        pc = p + 1;
    }
    B = InitialPositions.at(p);


    double WFa, WFb, WFc, WFd, WFe;// Values of Wave function in new points A/B/C
    WFa = WFValues.at(pa);
    WFb = WFValues.at(pb);
    WFc = WFValues.at(pc);
    D = (A+B)*0.5;
    E = (B+C)*0.5;
    PosX.at(0) = D;
    WFd = fabs ( WaveFunction(n, c, PosX, K, norm) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction(n, c, PosX, K, norm) );


for (int y = 0; y < NumOfBisections; y++)
{

//--------------- "D" ---------------
//    if ( WFd <= WFa ) // means f(D) <= f(A)
  //  {
        if ( WFd >= WFb  ) // means f(D) <= f(A), f(D) >= f(B), f(C) > f(B)
        {// new bracket (A=D, B, C)
            WFa = WFd;
            A = D;
        }
        else// means f(D) <= f(A), f(D) < f(B), f(C) > f(B)
        {// new bracket (A, B=D, C=B)
            WFc = WFb;
            WFb = WFd;
            C = B;
            B = D;
        }
    //}

//--------------- "E" ---------------
    //if ( WFe <= WFc ) // means f(E) <= f(C)
    //{
        if ( WFe >= WFb  ) // means f(E) <= f(C), f(E) >= f(B), f(A) > f(B)
        {// new bracket (A, B, C=E)
            WFc = WFe;
            C = E;
        }
        else// means f(E) <= f(C), f(E) < f(B), f(C) > f(B)
        {// new bracket (A=B, B=E, C)
            WFa = WFb;
            WFb = WFe;
            A = B;
            B = E;
        }
    //}
    D = (A+B)*0.5;
    E = (B+C)*0.5;
    PosX.at(0) = D;
    WFd = fabs ( WaveFunction(n, c, PosX, K, norm) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction(n, c, PosX, K, norm) );
}
  //  cout << "\n" << B << "   " << WFb << "\n";
    Minima.push_back( B );

}








void permutator::MinimaFinder_strong(vector<double> &Minima,
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
                              double norm)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The MinimaFinder function find the value of proper minimum of prob density based on the phase jump and value of prob
density in the point corresponding to the phase discontinuity. We know that the jump phase shouldn't be bigger than Pi :)

    Finding the values of phase and prob density in the InitialDivision number of points we can find an interesting
bracket. After that we use the bisection method. Phase and Prob density plot points obtained during initial calculation
may be used to prepare plots - so we can kill two birds with one stone.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

// ======================================================================================================================
// 1) Calculation of the Phase and Prob density values in initial points
// ======================================================================================================================

    vector<double>  PosX = Chain.at(ChainElement); // Choosing the collection of positions from the Chain (  .at(ChainElement)  )
    vector<double>  InitialPositions, PhaseValues, WFValues, SinValues, CosValues, temp;
    double  Phase, WF, SinValue, CosValue;

//Preparation of initial set of points (0,L*1/(InitialDiv-1),...,L)

    for(int i = 0; i < InitialDivision; i++)
    {
        InitialPositions.push_back( i*L/(InitialDivision - 1) );
    }

//Values of Phase and fabs(Wave Function)^2 in all the points from initial set

    for(int i = 0; i < InitialDivision; i++)
    {
        PosX.at(0) = InitialPositions.at(i);
        temp = PhaseCalculation(WaveFunction_strong(n, c, PosX, K, norm) );

        WF = temp.at(0);
        Phase = temp.at(1);
        SinValue = temp.at(2);
        CosValue = temp.at(3);

        WFValues.push_back( WF );        // Set of Prob denisty values in initial points
        PhaseValues.push_back( Phase );  // Set of Phase values in initial points
        SinValues.push_back( SinValue ); // Set of Sin values in initial points
        CosValues.push_back( CosValue ); // Set of Cos values in initial points

       // cout << i << "   " << PhaseValues.at(i) << "   " << WFValues.at(i) <<"\n";
    }

//Now we want to append values of Phase and Prob (and Sin and Cos) density to matrices

    AddCollection(PhaseMatrix, PhaseValues);
    AddCollection(WFMatrix, WFValues);

    AddCollection(SinMatrix, SinValues);
    AddCollection(CosMatrix, CosValues);

/*=====


// ======================================================================================================================
// 2) Proper value of the Phase Jump searching
// ======================================================================================================================

//Calculation of diferences between all 2 phase values: fabs( Phase(j+1)-Phase(j) ); and finding JUMPS!

    vector<double> JumpValue;
    vector<int> JumpPositions, JumpPositions2;
    double diff;

    for(int i = 0; i < InitialDivision-1; i++)
    {
        //cout << PhaseValues.at(i+1) << "   " << PhaseValues.at(i) <<  "   " <<   PhaseValues.at(i + 1) - PhaseValues.at( i )  <<"\n";
        diff =  fabs ( PhaseValues.at(i + 1) - PhaseValues.at( i ) );
        JumpValue.push_back( diff );

        //cout<< "\n diff: "<< diff << "\n";

        if( diff >= JumpRestriction)
        {
            JumpPositions.push_back(i);
        }
    }
/*
    cout << "\n" << "-----1------" << "\n";

    for(unsigned i = 0; i < JumpPositions.size(); i ++)
    {
    cout << JumpPositions.at(i) << "\n";
    }
*/

    /*=======


//We need to find the smallest value of Prob density in the set of points around the JUMPS

    int s, p;
    for(unsigned int j = 0; j < JumpPositions.size(); j++)
    {
        s = JumpPositions.at(j);
        JumpPositions2.push_back(s);
        JumpPositions2.push_back(s+1);
    }

/*
    cout << "-----2------" << "\n";
    for(unsigned i = 0; i < JumpPositions2.size(); i++)
    {
        cout << JumpPositions2.at(i) << "\n";
    }
*/



    /*=======


    // p is the number of element corresponding to minimal value
    // of wave function inside the bracket corresponding to the jump of phase

    p = JumpPositions2.at(0);
    for(int j = 1; j < JumpPositions2.size(); j++)
    {
        s = JumpPositions2.at(j);
        if( WFValues.at(p) < WFValues.at(s) )
        {
        }
        else{
            p = s;
        }
    }
    //cout << "\n" << "============" << "\n";
    //cout << p << "\n";


// ======================================================================================================================
// 3) Last step - BISECTION METHOD
// ======================================================================================================================

// (A=p-1,B=p,C=p+1) st. f(A)>f(B),f(C)>f(B) - bisection method

    double A,B,C,D,E;
    int pa, pb = p, pc; // number of collection element corresponding to position A/B/C

    if ( p > 0 )
    {
        A = InitialPositions.at(p-1);
        pa = p - 1;
        if ( p < InitialDivision-1 )
            {
                C = InitialPositions.at(p+1);
                pc = p + 1;
            }
            else //Means p = InitialDivision-1
            {
                C = InitialPositions.at(0);
                pc = 0;
            }
    }else//Means p = 0
    {
        A = InitialPositions.at(InitialDivision-1);
        pa = InitialDivision - 1;
        C = InitialPositions.at(p+1);
        pc = p + 1;
    }
    B = InitialPositions.at(p);


    double WFa, WFb, WFc, WFd, WFe;// Values of Wave function in new points A/B/C
    WFa = WFValues.at(pa);
    WFb = WFValues.at(pb);
    WFc = WFValues.at(pc);
    D = (A+B)*0.5;
    E = (B+C)*0.5;
    PosX.at(0) = D;
    WFd = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );



for (int y = 0; y < NumOfBisections; y++)
{

//--------------- "D" ---------------
//    if ( WFd <= WFa ) // means f(D) <= f(A)
  //  {
        if ( WFd >= WFb  ) // means f(D) <= f(A), f(D) >= f(B), f(C) > f(B)
        {// new bracket (A=D, B, C)
            WFa = WFd;
            A = D;
        }
        else// means f(D) <= f(A), f(D) < f(B), f(C) > f(B)
        {// new bracket (A, B=D, C=B)
            WFc = WFb;
            WFb = WFd;
            C = B;
            B = D;
        }
    //}

//--------------- "E" ---------------
    //if ( WFe <= WFc ) // means f(E) <= f(C)
    //{
        if ( WFe >= WFb  ) // means f(E) <= f(C), f(E) >= f(B), f(A) > f(B)
        {// new bracket (A, B, C=E)
            WFc = WFe;
            C = E;
        }
        else// means f(E) <= f(C), f(E) < f(B), f(C) > f(B)
        {// new bracket (A=B, B=E, C)
            WFa = WFb;
            WFb = WFe;
            A = B;
            B = E;
        }
    //}
    D = (A+B)*0.5;
    E = (B+C)*0.5;
    PosX.at(0) = D;
    WFd = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );
}


  //  cout << "\n" << B << "   " << WFb << "\n";
    Minima.push_back( B );



=====*/
}





void permutator::MinimaFinder_unnormalized(vector<double> &Minima,
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
                              double norm)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The MinimaFinder function find the value of proper minimum of prob density based on the phase jump and value of prob
density in the point corresponding to the phase discontinuity. We know that the jump phase shouldn't be bigger than Pi :)

    Finding the values of phase and prob density in the InitialDivision number of points we can find an interesting
bracket. After that we use the bisection method. Phase and Prob density plot points obtained during initial calculation
may be used to prepare plots - so we can kill two birds with one stone.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

// ======================================================================================================================
// 1) Calculation of the Phase and Prob density values in initial points
// ======================================================================================================================

    vector<double>  PosX = Chain.at(ChainElement); // Choosing the collection of positions from the Chain (  .at(ChainElement)  )
    vector<double>  InitialPositions, PhaseValues, WFValues, SinValues, CosValues, temp;
    double  Phase, WF, SinValue, CosValue;

//Preparation of initial set of points (0,L*1/(InitialDiv-1),...,L)

    for(int i = 0; i < InitialDivision; i++)
    {
        InitialPositions.push_back( i*L/(InitialDivision - 1) );
    }

//Values of Phase and fabs(Wave Function)^2 in all the points from initial set

    for(int i = 0; i < InitialDivision; i++)
    {
        PosX.at(0) = InitialPositions.at(i);
        temp = PhaseCalculation(WaveFunction_unnormalized(n, c, PosX, K, norm) );

        WF = temp.at(0);
        Phase = temp.at(1);
        SinValue = temp.at(2);
        CosValue = temp.at(3);

        WFValues.push_back( WF );        // Set of Prob denisty values in initial points
        PhaseValues.push_back( Phase );  // Set of Phase values in initial points
        SinValues.push_back( SinValue ); // Set of Sin values in initial points
        CosValues.push_back( CosValue ); // Set of Cos values in initial points

       // cout << i << "   " << PhaseValues.at(i) << "   " << WFValues.at(i) <<"\n";
    }

//Now we want to append values of Phase and Prob (and Sin and Cos) density to matrices

    AddCollection(PhaseMatrix, PhaseValues);
    AddCollection(WFMatrix, WFValues);

    AddCollection(SinMatrix, SinValues);
    AddCollection(CosMatrix, CosValues);

/*=====


// ======================================================================================================================
// 2) Proper value of the Phase Jump searching
// ======================================================================================================================

//Calculation of diferences between all 2 phase values: fabs( Phase(j+1)-Phase(j) ); and finding JUMPS!

    vector<double> JumpValue;
    vector<int> JumpPositions, JumpPositions2;
    double diff;

    for(int i = 0; i < InitialDivision-1; i++)
    {
        //cout << PhaseValues.at(i+1) << "   " << PhaseValues.at(i) <<  "   " <<   PhaseValues.at(i + 1) - PhaseValues.at( i )  <<"\n";
        diff =  fabs ( PhaseValues.at(i + 1) - PhaseValues.at( i ) );
        JumpValue.push_back( diff );

        //cout<< "\n diff: "<< diff << "\n";

        if( diff >= JumpRestriction)
        {
            JumpPositions.push_back(i);
        }
    }
/*
    cout << "\n" << "-----1------" << "\n";

    for(unsigned i = 0; i < JumpPositions.size(); i ++)
    {
    cout << JumpPositions.at(i) << "\n";
    }
*/

    /*=======


//We need to find the smallest value of Prob density in the set of points around the JUMPS

    int s, p;
    for(unsigned int j = 0; j < JumpPositions.size(); j++)
    {
        s = JumpPositions.at(j);
        JumpPositions2.push_back(s);
        JumpPositions2.push_back(s+1);
    }

/*
    cout << "-----2------" << "\n";
    for(unsigned i = 0; i < JumpPositions2.size(); i++)
    {
        cout << JumpPositions2.at(i) << "\n";
    }
*/



    /*=======


    // p is the number of element corresponding to minimal value
    // of wave function inside the bracket corresponding to the jump of phase

    p = JumpPositions2.at(0);
    for(int j = 1; j < JumpPositions2.size(); j++)
    {
        s = JumpPositions2.at(j);
        if( WFValues.at(p) < WFValues.at(s) )
        {
        }
        else{
            p = s;
        }
    }
    //cout << "\n" << "============" << "\n";
    //cout << p << "\n";


// ======================================================================================================================
// 3) Last step - BISECTION METHOD
// ======================================================================================================================

// (A=p-1,B=p,C=p+1) st. f(A)>f(B),f(C)>f(B) - bisection method

    double A,B,C,D,E;
    int pa, pb = p, pc; // number of collection element corresponding to position A/B/C

    if ( p > 0 )
    {
        A = InitialPositions.at(p-1);
        pa = p - 1;
        if ( p < InitialDivision-1 )
            {
                C = InitialPositions.at(p+1);
                pc = p + 1;
            }
            else //Means p = InitialDivision-1
            {
                C = InitialPositions.at(0);
                pc = 0;
            }
    }else//Means p = 0
    {
        A = InitialPositions.at(InitialDivision-1);
        pa = InitialDivision - 1;
        C = InitialPositions.at(p+1);
        pc = p + 1;
    }
    B = InitialPositions.at(p);


    double WFa, WFb, WFc, WFd, WFe;// Values of Wave function in new points A/B/C
    WFa = WFValues.at(pa);
    WFb = WFValues.at(pb);
    WFc = WFValues.at(pc);
    D = (A+B)*0.5;
    E = (B+C)*0.5;
    PosX.at(0) = D;
    WFd = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );



for (int y = 0; y < NumOfBisections; y++)
{

//--------------- "D" ---------------
//    if ( WFd <= WFa ) // means f(D) <= f(A)
  //  {
        if ( WFd >= WFb  ) // means f(D) <= f(A), f(D) >= f(B), f(C) > f(B)
        {// new bracket (A=D, B, C)
            WFa = WFd;
            A = D;
        }
        else// means f(D) <= f(A), f(D) < f(B), f(C) > f(B)
        {// new bracket (A, B=D, C=B)
            WFc = WFb;
            WFb = WFd;
            C = B;
            B = D;
        }
    //}

//--------------- "E" ---------------
    //if ( WFe <= WFc ) // means f(E) <= f(C)
    //{
        if ( WFe >= WFb  ) // means f(E) <= f(C), f(E) >= f(B), f(A) > f(B)
        {// new bracket (A, B, C=E)
            WFc = WFe;
            C = E;
        }
        else// means f(E) <= f(C), f(E) < f(B), f(C) > f(B)
        {// new bracket (A=B, B=E, C)
            WFa = WFb;
            WFb = WFe;
            A = B;
            B = E;
        }
    //}
    D = (A+B)*0.5;
    E = (B+C)*0.5;
    PosX.at(0) = D;
    WFd = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction_strong(n, c, PosX, K, norm) );
}


  //  cout << "\n" << B << "   " << WFb << "\n";
    Minima.push_back( B );



=====*/
}





//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//***************************** SCALAR PRODUCT - NUMERICAL INTEGRATION ***************************************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------


void permutator::IntegrationCollections(vector<double> K1,
                                        vector<double> K2,
                                        vector<vector<double>> &AllX,
                                        vector<double> &AllWs,
                                        vector<double> &AllNorms,
                                        double L,
                                        int PolynomialDegree,
                                        int IntervalDivision)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Preparation of collections which are indispensable to numerical integration (scalar product calculation <f2|f1>
over M-cube, where M = | K1.size() - K2.size() |.

Of course Measured_X.size() + M = N

In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int N1 = K1.size();
    int N2 = K2.size();
    int M  = N2 - abs(N1 - N2); // number of variables dx_1 dx_2 ... dx_M

// ======================================================================================================================
// 1) Integration intervals (+ coefficients) preparation
// ======================================================================================================================

    double b0 = L;
    vector<double> a;
    vector<double> b;
    vector<double> A1;
    vector<double> A2;
    double tmpA, tmpB;

    for (int i = 0; i < IntervalDivision; i++)
    {
        tmpA = ( b0/IntervalDivision )*i;
        tmpB = ( b0/IntervalDivision )*(i+1);

        a.push_back( tmpA );
        b.push_back( tmpB );

        A1.push_back( (tmpB - tmpA)/2. );
        A2.push_back( (tmpB + tmpA)/2. );
    }


// ======================================================================================================================
// 2) Weights and roots of the Legendre polynomial for PolynomialDegree
// We consider "only" PolynomialDegree in [5, 15]
// ======================================================================================================================

    vector<vector<double>> X_Roots = { {-0.90618, -0.538469, 0., 0.538469, 0.90618},       {-0.93247, -0.661209, -0.238619, 0.238619, 0.661209, 0.93247}, {-0.949108, -0.741531, -0.405845, 0., 0.405845, 0.741531, 0.949108}, {-0.96029, -0.796666, -0.525532, -0.183435, 0.183435, 0.525532, 0.796666, 0.96029}, {-0.96816, -0.836031, -0.613371, -0.324253, 0., 0.324253, 0.613371, 0.836031, 0.96816},      {-0.973907, -0.865063, -0.67941, -0.433395, -0.148874, 0.148874, 0.433395, 0.67941, 0.865063, 0.973907}, {-0.978229, -0.887063, -0.730152, -0.519096, -0.269543, 0., 0.269543, 0.519096, 0.730152, 0.887063, 0.978229}, {-0.981561, -0.904117, -0.769903, -0.587318, -0.367831, -0.125233, 0.125233, 0.367831, 0.587318, 0.769903, 0.904117, 0.981561}, {-0.984183, -0.917598, -0.801578, -0.642349, -0.448493, -0.230458, 0., 0.230458, 0.448493, 0.642349, 0.801578, 0.917598, 0.984183},   {-0.986284, -0.928435, -0.827201, -0.687293, -0.515249, -0.319112, -0.108055, 0.108055, 0.319112, 0.515249, 0.687293, 0.827201, 0.928435, 0.986284}, {-0.987993, -0.937273, -0.848207, -0.724418, -0.570972, -0.394151, -0.201194, 0., 0.201194, 0.394151, 0.570972, 0.724418, 0.848207, 0.937273, 0.987993} };

    vector<vector<double>> Weights = { {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}, {0.171324, 0.360762, 0.467914, 0.467914, 0.360762, 0.171324},  {0.129485, 0.279705, 0.38183, 0.417959, 0.38183, 0.279705, 0.129485},{0.101229, 0.222381, 0.313707, 0.362684, 0.362684, 0.313707, 0.222381, 0.101229},   {0.0812744, 0.180648, 0.260611, 0.312347, 0.330239, 0.312347, 0.260611, 0.180648, 0.0812744},{0.0666713, 0.149451, 0.219086, 0.269267, 0.295524, 0.295524, 0.269267, 0.219086, 0.149451, 0.0666713},  {0.0556686, 0.12558, 0.18629, 0.233194, 0.262805, 0.272925, 0.262805, 0.233194, 0.18629, 0.12558, 0.0556686},  {0.0471753, 0.106939, 0.160078, 0.203167, 0.233493, 0.249147, 0.249147, 0.233493, 0.203167, 0.160078, 0.106939, 0.0471753},     {0.040484, 0.0921215, 0.138874, 0.178146, 0.207816, 0.226283, 0.232552, 0.226283, 0.207816, 0.178146, 0.138874, 0.0921215, 0.040484}, {0.0351195, 0.0801581, 0.121519, 0.157203, 0.185538, 0.205198, 0.215264, 0.215264, 0.205198, 0.185538, 0.157203, 0.121519, 0.0801581, 0.0351195},    {0.0307532, 0.070366, 0.107159, 0.139571, 0.166269, 0.186161, 0.198431, 0.202578, 0.198431, 0.186161, 0.166269, 0.139571, 0.107159, 0.070366, 0.0307532} };

    vector<double> Xs = X_Roots.at(PolynomialDegree - 5);
    vector<double> Ws = Weights.at(PolynomialDegree - 5);


// ======================================================================================================================
// 3) Collection of positions + products of weights (all we need to perform summation)
// ======================================================================================================================

    vector<vector<double>> AllW;
    vector<vector<double>> AllNorm;
    vector<double> tmpAllX, tmpAllW ,tmpAllNorm;
    double tmpX;

    for(unsigned i1 = 0; i1 < A1.size(); i1++)
    {
        for(unsigned i2 = 0; i2 < Xs.size(); i2++)
        {
            tmpX = A1.at(i1)*Xs.at(i2) + A2.at(i1);
            tmpAllW.push_back( Ws.at(i2) );
            tmpAllNorm.push_back( A1.at(i1) );
            tmpAllX.push_back(tmpX);
        }
    }

    vector<vector<double>> CombX, CombW, CombNorm;
    vector<double> tmpCombX, tmpCombW, tmpCombNorm;
    double startX, startW, startNorm, CombSize1, CombSize2;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Preparation of all possible combinations
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(unsigned j1 = 0; j1 < tmpAllX.size(); j1++)
    {
        startX = tmpAllX.at(j1);
        startW = tmpAllW.at(j1);
        startNorm = tmpAllNorm.at(j1);

        tmpCombX.push_back(startX);
        tmpCombW.push_back(startW);
        tmpCombNorm.push_back(startNorm);

        AddCollection(CombX,tmpCombX);
        AddCollection(CombW,tmpCombW);
        AddCollection(CombNorm,tmpCombNorm);

        tmpCombX.clear();
        tmpCombW.clear();
        tmpCombNorm.clear();

//--------------------------------------------------------------------------------------------------------------------
        for(int m = 1; m < M; m++)
        {

            if( m == 1 )
            {
                CombSize1 = 0.;
                CombSize2 = 1.;
            }
            else
            {
                CombSize1 = CombSize2;
                CombSize2 += pow ( tmpAllX.size(), - 1. + m ) ;
            }

            for(int y = CombSize1; y < CombSize2; y++)
            {
                for(unsigned j2 = 0; j2 < tmpAllX.size(); j2++)
                {
                    tmpCombX = CombX.at(y);
                    tmpCombW = CombW.at(y);
                    tmpCombNorm = CombNorm.at(y);

                    //cout << "\n\n" << CombX[0][0];

                    tmpCombX.push_back(tmpAllX.at(j2));
                    tmpCombW.push_back(tmpAllW.at(j2));
                    tmpCombNorm.push_back(tmpAllNorm.at(j2));

                    AddCollection(CombX,tmpCombX);
                    AddCollection(CombW,tmpCombW);
                    AddCollection(CombNorm,tmpCombNorm);

                    tmpCombX.clear();
                    tmpCombW.clear();
                    tmpCombNorm.clear();
                }
            }
        }
//---------------------------------------------------------------------------------------------------------------------
            for(unsigned s = CombSize2; s < CombX.size(); s++ )
            {
                AddCollection(AllX, CombX.at(s));
                AddCollection(AllW, CombW.at(s));
                AddCollection(AllNorm, CombNorm.at(s));
            }

            CombX.clear();
            CombW.clear();
            CombNorm.clear();
    /*
            for(int u1 = 0; u1 < AllX.size(); u1++)
            {
                for(int u2 = 0; u2 < (AllX.at(u1)).size(); u2++)
                {
                    cout<< " , " << (AllX.at(u1)).at(u2);
                }
                cout<< "\n";
            }
    */

    }

        cout<<  "Number of relevant positions in "<< M<<"-cube:" << AllX.size() <<"\n";
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        double prod1, prod2;

        for(unsigned v1 = 0; v1 < AllW.size(); v1++)
        {
            prod1 = 1;
            prod2 = 1;
            for(unsigned v2 = 0; v2 < (AllW.at(v1)).size(); v2++)
            {
                prod1 = prod1 * (AllW.at(v1)).at(v2);
                prod2 = prod2 * (AllNorm.at(v1)).at(v2);
            }
            AllWs.push_back(prod1);
            AllNorms.push_back(prod2);
            //cout<< AllWs.at(v1) << " , " << AllNorms.at(v1) << "\n";
        }

}








void permutator::ScalarProduct(vector<double> K1,
                               vector<double> K2,
                               vector<vector<double>> AllX,
                               vector<double> AllWs,
                               vector<double> AllNorms,
                               double c,
                               vector<double> Measured_X,
                               vector<double> &ScalarProducts_re,
                               vector<double> &ScalarProducts_im,
                               double norm1,
                               double norm2)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Numerical integration (scalar product calculation <f2|f1> over M-cube, where M = | K1.size() - K2.size() |.
In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method - collections corresponding
to this method may be obtained using the function IntegrationCollections(...).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int N1 = K1.size();
    int N2 = K2.size();

// ======================================================================================================================
// 1) Last step of collections of positions preparation (here we assume that N2 > N1)
// ======================================================================================================================

    vector<vector<double>> X1, X2;
    vector<double> tmp;

    X1 = AllX;

    for(unsigned i = 0; i < AllX.size(); i++)
    {
        tmp = AllX.at(i);

        for(unsigned j = 0; j < Measured_X.size(); j++)
        {
            tmp.push_back( Measured_X.at(j) );
        }
        AddCollection(X2,tmp);
        tmp.clear();
    }

// ======================================================================================================================
// 2) Summation
// ======================================================================================================================

    Complex f0, f2;


    int s;
    double re,im;

    Complex int_value(0.,0.);
    Complex val2(0.,0.);

    for(s = 0; s < X1.size(); s++)
    {
        f0 = WaveFunction(N1, c, X1.at(s), K1, norm1);
        Complex f1( f0.real(), - f0.imag() );
        f2 = WaveFunction(N2, c, X2.at(s), K2, norm2);

        int_value = int_value + AllNorms.at(s) * AllWs.at(s) * f1 * f2;

        val2 = val2 +  AllNorms.at(s) * AllWs.at(s) * f0 * f1;


        if( s % 5000 ==0)
        {
            //cout<< s << " ";
        }
    }

    re=int_value.real();
    im=int_value.imag();

    Complex IntValue(re,im);

    ScalarProducts_re.push_back( IntValue.real() );
    ScalarProducts_im.push_back( IntValue.imag() );
    //cout<< "\n";


    cout<< val2.real() <<" , " << val2.imag() << "\n";

    //cout<< "(" << IntValue.real() << " , " << IntValue.imag() << ")" << "\n";


}






void permutator::RandPos_MonteCarloIntegration(int L,
                                               int NumOfParticles,
                                               int NumOfPositions,
                                               vector<vector<double>> &RandPos_Matrix)
{// the function randomly generates the collection of positions in M-Cube which we need
//to have in order to perform numerical integration using the Monte Carlo method

    double RandPos;

    for(int i = 0; i < NumOfPositions; i++)
    {
        vector<double> TmpVec;
        for(int j = 0; j < NumOfParticles; j++)
        {
            RandPos =  GetRandom(0,L); // Zmienilem z (0,L)
            TmpVec.push_back( RandPos );
        }
        RandPos_Matrix.push_back( TmpVec );
    }
}





void permutator::MonteCarlo_ScalarProduct(int Type,
                                          int L,
                                          vector<double> K1,
                                          vector<double> K2,
                                          vector<vector<double>> RandPos_Matrix,
                                          double c,
                                          vector<double> Measured_X,
                                          vector<double> &ScalarProducts_re,
                                          vector<double> &ScalarProducts_im,
                                          double norm1,
                                          double norm2)
{// We need to choose a type of method:
//  Type 1: Very easy method! Value = L^d (f(set X_1)/n + f(set X_2)/n ... f(set X_n)/n)
//  Type 2: Method based on the number of points counted "under the curve"

    int NumberOfPoints = RandPos_Matrix.size();
    int NumOfVariables = (RandPos_Matrix.at(0)).size();
    int N1 = K1.size();
    int N2 = K2.size();
    //here we assume that N2 > N1

    vector<vector<double>> X1, X2;
    vector<double> tmp;

    X1 = RandPos_Matrix;

    for(int i = 0; i < NumberOfPoints; i++)
    {
        tmp = RandPos_Matrix.at(i);

        for(unsigned j = 0; j < Measured_X.size(); j++)
        {
            tmp.push_back( Measured_X.at(j) );
        }

        AddCollection(X2,tmp);
        tmp.clear();
    }


    double re, im;

    if( Type == 1 )
    {

        Complex f0, f2;
        Complex Coefficient( pow(L, NumOfVariables) * 1./NumberOfPoints, 0. );
        Complex int_value(0.,0.);

        for(int s = 0; s < X1.size(); s++)
        {
            f0 = WaveFunction(N1, c, X1.at(s), K1, norm1);
            Complex f1( f0.real(), - f0.imag() );
            f2 = WaveFunction(N2, c, X2.at(s), K2, norm2);

            int_value = int_value +  Coefficient * f1 * f2;

            if( s % 5000 ==0)
            {
                //cout<< s << " ";
            }
        }

        re = int_value.real();
        im = int_value.imag();


    }
    else{

    Complex f0, f2, value;
    double ValRe, ValIm, RandRe, RandIm;
    double range = 3.;
    int NumRe = 0;
    int NumIm = 0;

        for(int s = 0; s < X1.size(); s++)
        {
            f0 = WaveFunction(N1, c, X1.at(s), K1, norm1);
            Complex f1( f0.real(), - f0.imag() );
            f2 = WaveFunction(N2, c, X2.at(s), K2, norm2);

            value = f1 * f2;
            ValRe = value.real();
            ValIm = value.imag();

            RandRe = GetRandom( -range , range);
            RandIm = GetRandom( -range , range);

//---------------- Re part -------------------------------------------------
            if( ValRe > 0. )
            {
                if( RandRe > 0. )
                {
                    if( ValRe > RandRe )
                    {
                        NumRe = NumRe + 1;
                    }
                }
            }
            else
            {
                if( RandRe < 0. )
                {
                    if( ValRe < RandRe )
                    {
                        NumRe = NumRe - 1;
                    }
                }
            }

//---------------- Im part -------------------------------------------------
            if( ValIm > 0. )
            {
                if( RandIm > 0. )
                {
                    if( ValIm > RandIm )
                    {
                        NumIm = NumIm + 1;
                    }
                }
            }
            else
            {
                if( RandIm < 0. )
                {
                    if( ValIm < RandIm )
                    {
                        NumIm = NumIm - 1;
                    }
                }
            }


        }

        re = 2. * range * pow( L ,NumOfVariables ) * (1.* NumRe)/NumberOfPoints;
        im = 2. * range * pow( L ,NumOfVariables ) * (1.* NumIm)/NumberOfPoints;

    }
    ScalarProducts_re.push_back( re );
    ScalarProducts_im.push_back( im );

}








void permutator::Scalar_GaussLegendre_Loop(int NumOfQuasimomentaSets,
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
                                           std::ofstream &sc_im)
{
//============== Scalar products calculation (Gauss-Legendre)  ====================================================================

    vector<vector<double>> AllX;
    vector<double> AllWs, AllNorms;


    IntegrationCollections(K_Matrix.at(0),K2,AllX,AllWs,AllNorms,L,
                           PolynomialDegree,IntervalDivision);


    vector<double>  ScalarProducts_re, ScalarProducts_im, Sc_re_tmp, Sc_im_tmp, K1;
    vector<int> Parallel_Order;
    Sc_re_tmp = vector<double> (NumOfQuasimomentaSets);
    Sc_im_tmp = vector<double> (NumOfQuasimomentaSets);


    int u = 0;
    int Par_i1; // iterator in the first parallelized loop
    double Par_norm2 = norm; // norm of the initial state
    double Par_norm1; // norm of tha Par_i1-th state of states considered after measurement


    #pragma omp parallel private(Par_i1 ,K1 ,Par_norm1)
    #pragma omp for
    for(Par_i1 = 0; Par_i1 < NumOfQuasimomentaSets; Par_i1++)
    {
        K1 = K_Matrix.at( Par_i1 );
        Par_norm1 = Norms.at( Par_i1 );
        ScalarProduct(K1, K2, AllX, AllWs, AllNorms,
                      c, Measured_X, ScalarProducts_re, ScalarProducts_im, Par_norm1, Par_norm2);

        Parallel_Order.push_back( Par_i1 );
        u++;
        if( u % 1 == 0)
        {
            cout<< u << " ";
        }
        //cout << "\n \n Row number of K_Matrix: " << Par_i1 << ", value of re part: " << ScalarProducts_re.back() << ", par order: " << Parallel_Order.back() << "\n \n";
    }

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
    sc_re.close();
    sc_im.close();



}






void permutator::Scalar_MonteCarlo_Loop(int Type,
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
                                        std::ofstream &sc_im)
{
    //============== Scalar products calculation (monte carlo)  ====================================================================


    vector<vector<double>> RandPos_Matrix;

    RandPos_MonteCarloIntegration(L, NumOfParticles, NumOfPositions, RandPos_Matrix);


    vector<double>  ScalarProducts_re, ScalarProducts_im, Sc_re_tmp, Sc_im_tmp, K1;
    vector<int> Parallel_Order;
    Sc_re_tmp = vector<double> (NumOfQuasimomentaSets);
    Sc_im_tmp = vector<double> (NumOfQuasimomentaSets);

    int u = 0;
    int Par_i1; // iterator in the first parallelized loop
    double Par_norm2 = norm; // norm of the initial state
    double Par_norm1; // norm of tha Par_i1-th state of states considered after measurement

    #pragma omp parallel private(Par_i1, K1 ,Par_norm1)
    #pragma omp for
    for(Par_i1 = 0; Par_i1 < NumOfQuasimomentaSets; Par_i1++)
    {
        K1 = K_Matrix.at( Par_i1 );
        Par_norm1 = Norms.at( Par_i1 );

        MonteCarlo_ScalarProduct(Type,L, K1, K2, RandPos_Matrix, c,
                                         Measured_X, ScalarProducts_re, ScalarProducts_im,
                                         Par_norm1, Par_norm2);

        Parallel_Order.push_back( Par_i1 );
        u++;
        if( u % 1 == 0)
        {
         cout << u << " " <<  "( " << ScalarProducts_re.back() << " , " << ScalarProducts_im.back() << " )" << "\n\n" ;
        }
    //cout << "\n \n Row number of K_Matrix: " << Par_i1 << ", value of re part: " << ScalarProducts_re.back() << ", par order: " << Parallel_Order.back() << "\n \n";
    }

    // ======== Sort and writing to file :) =================================================================================
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

    sc_re.close();
    sc_im.close();

}








//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//***************** RELEVANT COEFFICIENTS AND WAVE FUNCTIONS AFTER MEASUREMENT *******************************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

void permutator::RelevantCoefficients(vector<double> &ScalarProducts_re,
                                      vector<double> &ScalarProducts_im,
                                      vector<vector<double>> &K_Matrix,
                                      double CutOff_Coefficient,
                                      vector<double> &Norms)
{
    //We need to choose only important coefficients (in fact values of scalar products calculated before)
    //CutOff_Coefficient - the relevant coefficients are not smaller than MAX[abs(coefficients)]*CutOff_Coefficient

    int n = ScalarProducts_re.size();
    int PositionOfMax = 0;
    double tmp, MaxValue;
    vector<double> AbsValues, Tmp_re, Tmp_im, Tmp_K, Tmp_norm;
    vector<vector<double>> Tmp_K_Matrix;

    for(int i = 0; i < n; i++)
    {
        tmp = pow(   pow( ScalarProducts_re.at(i) , 2.) + pow( ScalarProducts_im.at(i) , 2.)   ,   0.5 );

        AbsValues.push_back( tmp );

        if( i == 0 )
        {}
        else
        {
            if( AbsValues.at( PositionOfMax ) >= tmp )
            {}
            else
            {
                PositionOfMax = i;
            }
        }
    }


    MaxValue = AbsValues.at( PositionOfMax );
    double NotLessThan = MaxValue * CutOff_Coefficient; // "real" cut off

    for(int j = 0; j < n; j++)
    {
        if( AbsValues.at(j) >= NotLessThan )
        {
            Tmp_re.push_back( ScalarProducts_re.at(j) );
            Tmp_im.push_back( ScalarProducts_im.at(j) );
            Tmp_norm.push_back( Norms.at(j) );

            Tmp_K = K_Matrix.at(j);
            Tmp_K_Matrix.push_back( Tmp_K );


        }
    }

    // NEW RELEVANT VALUES OF SCALAR PRODUCTS
    ScalarProducts_re = Tmp_re;
    ScalarProducts_im = Tmp_im;

    //AND CORRESPONDING K_Matrix
    K_Matrix = Tmp_K_Matrix;
    Norms = Tmp_norm;
}







std::complex<double> permutator::WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix,
                                                              vector<double> ScalarProducts_re,
                                                              vector<double> ScalarProducts_im,
                                                              vector<double> Remaining_X,
                                                              double c,
                                                              double t,
                                                              vector<double> Norms)
{
// The wave function after measurement of particle positions - it is given by a sum over N-M particle WaveFunctions with proper coefficients (scalar products)
// Cause of the structure of the Wave Function we obtain time evolution in very easy way.

    int n = Remaining_X.size();
    int NumOfRows = K_Matrix.size();

    double re = 0., im = 0.;
    int i;

//#pragma omp parallel  private(i) reduction(+:re) reduction(+:im)
//{
    Complex WF(0. , 0.);
  //  #pragma omp for

    for(i = 0; i < NumOfRows; i++)
    {
        double energy = 0;

        for(int j = 0; j < n; j++)
        {
            energy = energy + pow( (K_Matrix.at(i)).at(j) , 2.);
        }

        Complex Coefficient( ScalarProducts_re.at(i) , ScalarProducts_im.at(i) );
        Complex I_Exp_t(0. , energy * t);
        Complex EvolutionExp;
        EvolutionExp = exp(- I_Exp_t);

        WF = WF +  EvolutionExp * Coefficient  * WaveFunction(n, c, Remaining_X, K_Matrix.at( i ), Norms.at(i) );

        //cout<< "\n time t=" << t << "\n\n";
        //cout<< "\n " << fabs( WF ) <<"    (" << EvolutionExp.real() << " , " <<  EvolutionExp.imag() <<  ")" << "\n\n";

    }

    re=WF.real();
    im=WF.imag();
//}

    Complex Wavefunction(re,im);

    return Wavefunction;

}






void permutator::Metropolis_TimeDependent(int steps,
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
                                          vector<double> Norms)
{
    double ratio;
    double dif;
    double random0to1;
    double wf0;
    double wf1;

    wf0 =  pow( fabs( WaveFunctionAfterMeasurement(K_Matrix, ScalarProducts_re, ScalarProducts_im, X0, c, t, Norms)  ), 2.);

for(int u = 0; u < steps; u++)
{

    for(int i = 0; i < N; i++)
    {
        Positions(X0, X1, delta, i, L);
    }
        wf1 =  pow( fabs( WaveFunctionAfterMeasurement(K_Matrix, ScalarProducts_re, ScalarProducts_im, X1, c, t, Norms)  ), 2.);

        ratio = wf1/wf0;

        if( ratio >= 1 )
        {
            WaveFValues.push_back(wf1);
            AddCollection(Chain,X1);
        }
        else
        {
            random0to1 = GetRandom(0.,1.);
            dif = ratio - random0to1;
                if( dif > 0 )
                {
                    WaveFValues.push_back(wf1);
                    AddCollection(Chain,X0);
                }
                else
                {
                    X1 = X0;
                    wf1 = wf0;
                }
        }
     X0 = X1;
     wf0 = wf1;

     if (u % 100 == 0)
     {
         cout << u << ", ";
     }

}

}










//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//******************** FUNCTIONS FOR THE REPULSIVE CASE (COMPLEX QUASIMOMENTA) *******************************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
// ###################################################################################################################################
// R like Repulsive is added to the name of function in the following way R_FunctionName
// Now K -> K_re (realis of K) and K_im (imaginalis of K)
// ###################################################################################################################################





// ###################################################################################################################################
// =============================== R_sProduct ========================================================================================
// ###################################################################################################################################


std::complex<double> permutator::R_sProduct(int n,
                                            int PermNumber,
                                            double c,
                                            vector<double> X,
                                            vector<double> K_re,
                                            vector<double> K_im)
{
// ====================================================================================================================================
// 1) Product (j>k which appears in the wave function) definition. PermNumber is the number of permutation in permutations array (row)
// ====================================================================================================================================

    Complex  prod(1.,0.);
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < j; k++)
        {
        prod = prod* Complex( K_re[GetPermEl(PermNumber,j)] - K_re[GetPermEl(PermNumber,k)], K_im[GetPermEl(PermNumber,j)] - K_im[GetPermEl(PermNumber,k)] -c*SignFunction( X[j], X[k]) ) ;
        }
    }
    // cout << prod.real() << ',' << prod.imag();


// ====================================================================================================================================
// 2) Exponent with a sum over X[]K[]
// ====================================================================================================================================

    Complex sumexp (0.,0.);
    for(int j = 0; j < n; j++)
    {
        sumexp = sumexp + Complex ( K_re[GetPermEl(PermNumber,j)] * X[j] , K_im[GetPermEl(PermNumber,j)] * X[j]);
    }
    Complex ImSumExp(- sumexp.imag() , sumexp.real() );
    Complex ImExp;
    ImExp = exp(ImSumExp);
    // cout << ImExp.real() << ',' << ImExp.imag();

// ====================================================================================================================================
// 3) Sign of the permutation given by PermNumber
// ====================================================================================================================================

    int SPV = _PermSigns.at(PermNumber); // + or - 1 (Parity)
    Complex SignValue(SPV, 0.);

// ====================================================================================================================================
// 4) Finally, we can calculate one element of "big" sum apearing in the wave function
// ====================================================================================================================================

    Complex ProductElement;
    ProductElement = SignValue * ImExp * prod;

/*
    cout <<'('<<ProductElement.real()<<','<<ProductElement.imag()<<')' ;
    cout <<'('<<SignValue.real()<<','<<SignValue.imag()<<')' ;
    cout << prod.real() << ImExp.real();
    cout << prod.imag() << ImExp.imag();
*/

   // cout << ProductElement.real() << ',' << ProductElement.imag();
    return ProductElement;
}





// ###################################################################################################################################
// =============================== R_WaveFunction ====================================================================================
// ###################################################################################################################################

std::complex<double> permutator::R_WaveFunction(int n,
                                                double c,
                                                vector<double> X,
                                                vector<double> K_re,
                                                vector<double> K_im,
                                                double normalizacja)
{
    // N! calculation - we need to know how many permutations we will have
    //*************************************************************************
            int Nfact = 1; // Nfact = N! (Nfactorial)
            for (int i = 2; i <= n; i++) // Nfactorial value (N!) calculation
            { Nfact = Nfact * i; }
            //cout << Nfact;
    // Actually, it is the number of rows of "perms"
    // ************************************************************************


    double re = 0., im = 0.;
    int i;
    //#pragma omp parallel  private(i) reduction(+:re) reduction(+:im)
    //{
        Complex wavefunction(0.,0.);
        //#pragma omp for
        for(i = 0; i < Nfact; i++)
        {
            wavefunction = wavefunction + R_sProduct(n, i, c, X, K_re, K_im);
        }
        re = wavefunction.real();
        im = wavefunction.imag();
    //}
    Complex Wavefunction(re,im);


    // Normalization factor

    Complex norm (Nfact, 0.);
    for (int j = 0; j < n ; j++)
    {
        for (int k = 0; k < j; k++)
        {
        Complex Kdiff( (K_re[j] - K_re[k]), ( K_im[j] - K_im[k] ) );
        norm = norm * ( Kdiff * Kdiff + Complex( c*c, 0. ) );
        }
    }


    Complex InvNorm( norm.real()/( pow( norm.real(),2.) + pow( norm.imag(),2.) ) , - norm.imag()/( pow( norm.real(),2.) + pow( norm.imag(),2.) ) );

    Complex Norm;
    Norm = std::sqrt( InvNorm );
    Complex Normalizacja(1./sqrt(normalizacja),0.);
    Wavefunction = Normalizacja * Norm *  Wavefunction;

    //cout << fabs(wavefunction) <<  "\n";
    //cout << Norm.real() << ", " << Norm.imag() << "\n";
    //cout << Wavefunction.real() << ',' << Wavefunction.imag() << "\n";

    //wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return Wavefunction;
}







double permutator::R_TwoBodyWF_approx_merged(double x1, double x2,double L, double c)
// Repulsive wave function of ground state in case of 2 particles. I froced the ring topology!
// It is enought to force periodicity on this 2-particle wave function because full approximated
// wave function is 2-particle reducible.
{
    double value;

    if( x1 < L/2.)
    {
        if( x2 < x1 + L/2.)
        {
            value = exp(c/2. * fabs( x1 - x2 ) );
        }
        else
        {
            value = exp(c/2. * fabs( L + x1 - x2 ) );
        }
    }
    else
    {
        if( x2 > x1 - L/2.)
        {
            value = exp(c/2. * fabs( x1 - x2 ) );
        }
        else
        {
            value = exp(c/2. * fabs( ( L - x1 ) + x2 ) );
        }
    }

    return value;
    
}






double permutator::R_WF_approx(vector<double> X, double L, double c)
//Approximated wave function (see eq (85) Y. Castin, C. Herzog "Bose Einstein condensates in symmetry breaking states"
// - arXiv:cond-mat/0012040v2
// with forced periodicity :D
{
    int N = X.size();

// ======== calculation of ~normalization (see eq. (86)) ===============================

    int n = 1;

    for(int i = 1; i < N; i++)
    {
        n = n*i;
    }

    double norm = pow(  (1.* n)/(L * N), -1.) * pow( .5 * c/2, -1.* (N-1) ) ;

// ======== calculation of N-body wave function ========================================

    double value = 1.;

    for(int i = 0; i < N; i++)
    {
        for(int j = i + 1; j < N; j++)
        {
            value = value * R_TwoBodyWF_approx_merged(X[i], X[j], L, c);
        }
    }

    value = value * norm;

    return value /* pow( 2.93202 *pow(10. , - 13.) ,-.5)*/;

}






Complex permutator::R_pushed_WF_approx(vector<double> X, double P, double L, double c)
// Approximated wave function of ground state with additional momentum of centre of mass (P)
{
    int N = X.size();
    double x = 0.;
    for(int i = 0; i < N; i++)
    {
        x = x + X[i];
    }

    Complex ImPX(0. , P * x);
    Complex ExpP;

    ExpP = exp(ImPX);

    Complex value0(R_WF_approx(X, L, c), 0.);
    Complex value;

    return value = ExpP * value0  * Complex(pow(2.87115 * pow(10. , -9.), -.5) , 0. );
//*/
}







void permutator::R_Metropolis_approx(int OutSteps,
                                   int InsideSteps,
                                   int N,
                                   double c,
                                   vector<vector<double> > &Chain,
                                   vector<double> &WaveFValues,
                                   vector<double> &X1,
                                   double delta,
                                   double L
                                   )
{
    double ratio;
    double dif;
    double random0to1;
    double wf0;
    double wf1;

for(int f = 0; f < OutSteps; f++)
{
    vector<double> X0;
    for(int z = 0; z < N; z++)
    {
        X0.push_back(GetRandom(0.,1.));
    }

    wf0 =  pow (fabs(R_WF_approx(X0, L, c) ), 2.);

for(int u = 0; u < InsideSteps; u++)
{


    for(int i = 0; i < N; i++)
    {
        Positions(X0, X1, delta, i, L);
    }
        wf1 =  pow (fabs(R_WF_approx(X1, L, c) ), 2.);

        ratio = wf1/wf0;

        if( ratio >= 1. )
        {
            WaveFValues.push_back(wf1);
            AddCollection(Chain,X1);
            X0 = X1;
            wf0 = wf1;
        }
        else
        {
            random0to1 = GetRandom(0.,1.);
            dif = ratio - random0to1;
                if( dif > 0. )
                {
                    WaveFValues.push_back(wf1);
                    AddCollection(Chain,X1);
                    X0 = X1;
                    wf0 = wf1;
                }
                else
                {
                    X1 = X0;
                    wf1 = wf0;
                }
        }

     if(f%100 == 0)
     {
        cout<< "\r"<< f << " ";
     }
}

}


}










// ###################################################################################################################################
// =============================== R_IntegrationCollections ==========================================================================
// ###################################################################################################################################


void permutator::R_IntegrationCollections(vector<double> K1_re,
                                          vector<double> K2_re,
                                          vector<vector<double>> &AllX,
                                          vector<double> &AllWs,
                                          vector<double> &AllNorms,
                                          double L,
                                          int PolynomialDegree,
                                          int IntervalDivision)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Preparation of collections which are indispensable to numerical integration (scalar product calculation <f2|f1>
over M-cube, where M = | K1.size() - K2.size() |.

Of course Measured_X.size() + M = N

In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int N1 = K1_re.size();
    int N2 = K2_re.size();
    int M  = N2 - abs(N1 - N2); // number of variables dx_1 dx_2 ... dx_M

// ======================================================================================================================
// 1) Integration intervals (+ coefficients) preparation
// ======================================================================================================================

    double b0 = L;
    vector<double> a;
    vector<double> b;
    vector<double> A1;
    vector<double> A2;
    double tmpA, tmpB;

    for (int i = 0; i < IntervalDivision; i++)
    {
        tmpA = ( b0/IntervalDivision )*i;
        tmpB = ( b0/IntervalDivision )*(i+1);

        a.push_back( tmpA );
        b.push_back( tmpB );

        A1.push_back( (tmpB - tmpA)/2. );
        A2.push_back( (tmpB + tmpA)/2. );
    }


// ======================================================================================================================
// 2) Weights and roots of the Legendre polynomial for PolynomialDegree
// We consider "only" PolynomialDegree in [5, 15]
// ======================================================================================================================

    vector<vector<double>> X_Roots = { {-0.90618, -0.538469, 0., 0.538469, 0.90618},       {-0.93247, -0.661209, -0.238619, 0.238619, 0.661209, 0.93247}, {-0.949108, -0.741531, -0.405845, 0., 0.405845, 0.741531, 0.949108}, {-0.96029, -0.796666, -0.525532, -0.183435, 0.183435, 0.525532, 0.796666, 0.96029}, {-0.96816, -0.836031, -0.613371, -0.324253, 0., 0.324253, 0.613371, 0.836031, 0.96816},      {-0.973907, -0.865063, -0.67941, -0.433395, -0.148874, 0.148874, 0.433395, 0.67941, 0.865063, 0.973907}, {-0.978229, -0.887063, -0.730152, -0.519096, -0.269543, 0., 0.269543, 0.519096, 0.730152, 0.887063, 0.978229}, {-0.981561, -0.904117, -0.769903, -0.587318, -0.367831, -0.125233, 0.125233, 0.367831, 0.587318, 0.769903, 0.904117, 0.981561}, {-0.984183, -0.917598, -0.801578, -0.642349, -0.448493, -0.230458, 0., 0.230458, 0.448493, 0.642349, 0.801578, 0.917598, 0.984183},   {-0.986284, -0.928435, -0.827201, -0.687293, -0.515249, -0.319112, -0.108055, 0.108055, 0.319112, 0.515249, 0.687293, 0.827201, 0.928435, 0.986284}, {-0.987993, -0.937273, -0.848207, -0.724418, -0.570972, -0.394151, -0.201194, 0., 0.201194, 0.394151, 0.570972, 0.724418, 0.848207, 0.937273, 0.987993} };

    vector<vector<double>> Weights = { {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}, {0.171324, 0.360762, 0.467914, 0.467914, 0.360762, 0.171324},  {0.129485, 0.279705, 0.38183, 0.417959, 0.38183, 0.279705, 0.129485},{0.101229, 0.222381, 0.313707, 0.362684, 0.362684, 0.313707, 0.222381, 0.101229},   {0.0812744, 0.180648, 0.260611, 0.312347, 0.330239, 0.312347, 0.260611, 0.180648, 0.0812744},{0.0666713, 0.149451, 0.219086, 0.269267, 0.295524, 0.295524, 0.269267, 0.219086, 0.149451, 0.0666713},  {0.0556686, 0.12558, 0.18629, 0.233194, 0.262805, 0.272925, 0.262805, 0.233194, 0.18629, 0.12558, 0.0556686},  {0.0471753, 0.106939, 0.160078, 0.203167, 0.233493, 0.249147, 0.249147, 0.233493, 0.203167, 0.160078, 0.106939, 0.0471753},     {0.040484, 0.0921215, 0.138874, 0.178146, 0.207816, 0.226283, 0.232552, 0.226283, 0.207816, 0.178146, 0.138874, 0.0921215, 0.040484}, {0.0351195, 0.0801581, 0.121519, 0.157203, 0.185538, 0.205198, 0.215264, 0.215264, 0.205198, 0.185538, 0.157203, 0.121519, 0.0801581, 0.0351195},    {0.0307532, 0.070366, 0.107159, 0.139571, 0.166269, 0.186161, 0.198431, 0.202578, 0.198431, 0.186161, 0.166269, 0.139571, 0.107159, 0.070366, 0.0307532} };

    vector<double> Xs = X_Roots.at(PolynomialDegree - 5);
    vector<double> Ws = Weights.at(PolynomialDegree - 5);


// ======================================================================================================================
// 3) Collection of positions + products of weights (all we need to perform summation)
// ======================================================================================================================

    vector<vector<double>> AllW;
    vector<vector<double>> AllNorm;
    vector<double> tmpAllX, tmpAllW ,tmpAllNorm;
    double tmpX;

    for(unsigned i1 = 0; i1 < A1.size(); i1++)
    {
        for(unsigned i2 = 0; i2 < Xs.size(); i2++)
        {
            tmpX = A1.at(i1)*Xs.at(i2) + A2.at(i1);
            tmpAllW.push_back( Ws.at(i2) );
            tmpAllNorm.push_back( A1.at(i1) );
            tmpAllX.push_back(tmpX);
        }
    }

    vector<vector<double>> CombX, CombW, CombNorm;
    vector<double> tmpCombX, tmpCombW, tmpCombNorm;
    double startX, startW, startNorm, CombSize1, CombSize2;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Preparation of all possible combinations
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(unsigned j1 = 0; j1 < tmpAllX.size(); j1++)
    {
        startX = tmpAllX.at(j1);
        startW = tmpAllW.at(j1);
        startNorm = tmpAllNorm.at(j1);

        tmpCombX.push_back(startX);
        tmpCombW.push_back(startW);
        tmpCombNorm.push_back(startNorm);

        AddCollection(CombX,tmpCombX);
        AddCollection(CombW,tmpCombW);
        AddCollection(CombNorm,tmpCombNorm);

        tmpCombX.clear();
        tmpCombW.clear();
        tmpCombNorm.clear();

//------------------------------------------------------------------------------------------------------------------
        for(int m = 1; m < M; m++)
        {

            if( m == 1 )
            {
                CombSize1 = 0.;
                CombSize2 = 1.;
            }
            else
            {
                CombSize1 = CombSize2;
                CombSize2 += pow ( tmpAllX.size(), - 1. + m ) ;
            }

            for(int y = CombSize1; y < CombSize2; y++)
            {
                for(unsigned j2 = 0; j2 < tmpAllX.size(); j2++)
                {
                    tmpCombX = CombX.at(y);
                    tmpCombW = CombW.at(y);
                    tmpCombNorm = CombNorm.at(y);

                    tmpCombX.push_back(tmpAllX.at(j2));
                    tmpCombW.push_back(tmpAllW.at(j2));
                    tmpCombNorm.push_back(tmpAllNorm.at(j2));

                    AddCollection(CombX,tmpCombX);
                    AddCollection(CombW,tmpCombW);
                    AddCollection(CombNorm,tmpCombNorm);

                    tmpCombX.clear();
                    tmpCombW.clear();
                    tmpCombNorm.clear();
                }
            }
        }
//--------------------------------------------------------------------------------------------------------------------
            for(unsigned s = CombSize2; s < CombX.size(); s++ )
            {
                AddCollection(AllX, CombX.at(s));
                AddCollection(AllW, CombW.at(s));
                AddCollection(AllNorm, CombNorm.at(s));
            }

            CombX.clear();
            CombW.clear();
            CombNorm.clear();
    /*
            for(int u1 = 0; u1 < AllX.size(); u1++)
            {
                for(int u2 = 0; u2 < (AllX.at(u1)).size(); u2++)
                {
                    cout<< " , " << (AllX.at(u1)).at(u2);
                }
                cout<< "\n";
            }
    */

    }

        cout<<  "Number of relevant positions in M-cube:" << AllX.size() <<"\n";
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        double prod1, prod2;

        for(unsigned v1 = 0; v1 < AllW.size(); v1++)
        {
            prod1 = 1;
            prod2 = 1;
            for(unsigned v2 = 0; v2 < (AllW.at(v1)).size(); v2++)
            {
                prod1 = prod1 * (AllW.at(v1)).at(v2);
                prod2 = prod2 * (AllNorm.at(v1)).at(v2);
            }
            AllWs.push_back(prod1);
            AllNorms.push_back(prod2);
            //cout<< AllWs.at(v1) << " , " << AllNorms.at(v1) << "\n";
        }

}




void permutator::R_IntegrationCollections_approx(int N1,
                                                 int N2,
                                                 vector<vector<double>> &AllX,
                                                 vector<double> &AllWs,
                                                 vector<double> &AllNorms,
                                                 double L,
                                                 int PolynomialDegree,
                                                 int IntervalDivision)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Preparation of collections which are indispensable to numerical integration (scalar product calculation <f2|f1>
over M-cube, where M = | K1.size() - K2.size() |.

Of course Measured_X.size() + M = N

In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int M  = N2 - abs(N1 - N2); // number of variables dx_1 dx_2 ... dx_M

    cout << M << "\n\n";

// ======================================================================================================================
// 1) Integration intervals (+ coefficients) preparation
// ======================================================================================================================

    double b0 = L;
    vector<double> a;
    vector<double> b;
    vector<double> A1;
    vector<double> A2;
    double tmpA, tmpB;

    for (int i = 0; i < IntervalDivision; i++)
    {
        tmpA = ( b0/IntervalDivision )*i;
        tmpB = ( b0/IntervalDivision )*(i+1);

        a.push_back( tmpA );
        b.push_back( tmpB );

        A1.push_back( (tmpB - tmpA)/2. );
        A2.push_back( (tmpB + tmpA)/2. );
    }


// ======================================================================================================================
// 2) Weights and roots of the Legendre polynomial for PolynomialDegree
// We consider "only" PolynomialDegree in [5, 15]
// ======================================================================================================================

    vector<vector<double>> X_Roots = { {-0.90618, -0.538469, 0., 0.538469, 0.90618},       {-0.93247, -0.661209, -0.238619, 0.238619, 0.661209, 0.93247}, {-0.949108, -0.741531, -0.405845, 0., 0.405845, 0.741531, 0.949108}, {-0.96029, -0.796666, -0.525532, -0.183435, 0.183435, 0.525532, 0.796666, 0.96029}, {-0.96816, -0.836031, -0.613371, -0.324253, 0., 0.324253, 0.613371, 0.836031, 0.96816},      {-0.973907, -0.865063, -0.67941, -0.433395, -0.148874, 0.148874, 0.433395, 0.67941, 0.865063, 0.973907}, {-0.978229, -0.887063, -0.730152, -0.519096, -0.269543, 0., 0.269543, 0.519096, 0.730152, 0.887063, 0.978229}, {-0.981561, -0.904117, -0.769903, -0.587318, -0.367831, -0.125233, 0.125233, 0.367831, 0.587318, 0.769903, 0.904117, 0.981561}, {-0.984183, -0.917598, -0.801578, -0.642349, -0.448493, -0.230458, 0., 0.230458, 0.448493, 0.642349, 0.801578, 0.917598, 0.984183},   {-0.986284, -0.928435, -0.827201, -0.687293, -0.515249, -0.319112, -0.108055, 0.108055, 0.319112, 0.515249, 0.687293, 0.827201, 0.928435, 0.986284}, {-0.987993, -0.937273, -0.848207, -0.724418, -0.570972, -0.394151, -0.201194, 0., 0.201194, 0.394151, 0.570972, 0.724418, 0.848207, 0.937273, 0.987993} };

    vector<vector<double>> Weights = { {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}, {0.171324, 0.360762, 0.467914, 0.467914, 0.360762, 0.171324},  {0.129485, 0.279705, 0.38183, 0.417959, 0.38183, 0.279705, 0.129485},{0.101229, 0.222381, 0.313707, 0.362684, 0.362684, 0.313707, 0.222381, 0.101229},   {0.0812744, 0.180648, 0.260611, 0.312347, 0.330239, 0.312347, 0.260611, 0.180648, 0.0812744},{0.0666713, 0.149451, 0.219086, 0.269267, 0.295524, 0.295524, 0.269267, 0.219086, 0.149451, 0.0666713},  {0.0556686, 0.12558, 0.18629, 0.233194, 0.262805, 0.272925, 0.262805, 0.233194, 0.18629, 0.12558, 0.0556686},  {0.0471753, 0.106939, 0.160078, 0.203167, 0.233493, 0.249147, 0.249147, 0.233493, 0.203167, 0.160078, 0.106939, 0.0471753},     {0.040484, 0.0921215, 0.138874, 0.178146, 0.207816, 0.226283, 0.232552, 0.226283, 0.207816, 0.178146, 0.138874, 0.0921215, 0.040484}, {0.0351195, 0.0801581, 0.121519, 0.157203, 0.185538, 0.205198, 0.215264, 0.215264, 0.205198, 0.185538, 0.157203, 0.121519, 0.0801581, 0.0351195},    {0.0307532, 0.070366, 0.107159, 0.139571, 0.166269, 0.186161, 0.198431, 0.202578, 0.198431, 0.186161, 0.166269, 0.139571, 0.107159, 0.070366, 0.0307532} };

    vector<double> Xs = X_Roots.at(PolynomialDegree - 5);
    vector<double> Ws = Weights.at(PolynomialDegree - 5);


// ======================================================================================================================
// 3) Collection of positions + products of weights (all we need to perform summation)
// ======================================================================================================================

    vector<vector<double>> AllW;
    vector<vector<double>> AllNorm;
    vector<double> tmpAllX, tmpAllW ,tmpAllNorm;
    double tmpX;

    for(unsigned i1 = 0; i1 < A1.size(); i1++)
    {
        for(unsigned i2 = 0; i2 < Xs.size(); i2++)
        {
            tmpX = A1.at(i1)*Xs.at(i2) + A2.at(i1);
            tmpAllW.push_back( Ws.at(i2) );
            tmpAllNorm.push_back( A1.at(i1) );
            tmpAllX.push_back(tmpX);
        }
    }

    vector<vector<double>> CombX, CombW, CombNorm;
    vector<double> tmpCombX, tmpCombW, tmpCombNorm;
    double startX, startW, startNorm, CombSize1, CombSize2;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Preparation of all possible combinations
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(unsigned j1 = 0; j1 < tmpAllX.size(); j1++)
    {
        startX = tmpAllX.at(j1);
        startW = tmpAllW.at(j1);
        startNorm = tmpAllNorm.at(j1);

        tmpCombX.push_back(startX);
        tmpCombW.push_back(startW);
        tmpCombNorm.push_back(startNorm);

        AddCollection(CombX,tmpCombX);
        AddCollection(CombW,tmpCombW);
        AddCollection(CombNorm,tmpCombNorm);

        tmpCombX.clear();
        tmpCombW.clear();
        tmpCombNorm.clear();

//------------------------------------------------------------------------------------------------------------------
        for(int m = 1; m < M; m++)
        {

            if( m == 1 )
            {
                CombSize1 = 0.;
                CombSize2 = 1.;
            }
            else
            {
                CombSize1 = CombSize2;
                CombSize2 += pow ( tmpAllX.size(), - 1. + m ) ;
            }

            for(int y = CombSize1; y < CombSize2; y++)
            {
                for(unsigned j2 = 0; j2 < tmpAllX.size(); j2++)
                {
                    tmpCombX = CombX.at(y);
                    tmpCombW = CombW.at(y);
                    tmpCombNorm = CombNorm.at(y);

                    tmpCombX.push_back(tmpAllX.at(j2));
                    tmpCombW.push_back(tmpAllW.at(j2));
                    tmpCombNorm.push_back(tmpAllNorm.at(j2));

                    AddCollection(CombX,tmpCombX);
                    AddCollection(CombW,tmpCombW);
                    AddCollection(CombNorm,tmpCombNorm);

                    tmpCombX.clear();
                    tmpCombW.clear();
                    tmpCombNorm.clear();
                }
            }
        }
//--------------------------------------------------------------------------------------------------------------------
            for(unsigned s = CombSize2; s < CombX.size(); s++ )
            {
                AddCollection(AllX, CombX.at(s));
                AddCollection(AllW, CombW.at(s));
                AddCollection(AllNorm, CombNorm.at(s));
            }

            CombX.clear();
            CombW.clear();
            CombNorm.clear();
    /*
            for(int u1 = 0; u1 < AllX.size(); u1++)
            {
                for(int u2 = 0; u2 < (AllX.at(u1)).size(); u2++)
                {
                    cout<< " , " << (AllX.at(u1)).at(u2);
                }
                cout<< "\n";
            }
    */

    }

        cout<<  "Number of relevant positions in M-cube:" << AllX.size() <<"\n";
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        double prod1, prod2;

        for(unsigned v1 = 0; v1 < AllW.size(); v1++)
        {
            prod1 = 1;
            prod2 = 1;
            for(unsigned v2 = 0; v2 < (AllW.at(v1)).size(); v2++)
            {
                prod1 = prod1 * (AllW.at(v1)).at(v2);
                prod2 = prod2 * (AllNorm.at(v1)).at(v2);
            }
            AllWs.push_back(prod1);
            AllNorms.push_back(prod2);
            //cout<< AllWs.at(v1) << " , " << AllNorms.at(v1) << "\n";
        }

}













// ###################################################################################################################################
// =============================== R_ScalarProduct ===================================================================================
// ###################################################################################################################################


void permutator::R_ScalarProduct(vector<double> K1_re,
                                 vector<double> K1_im,
                                 vector<double> K2_re,
                                 vector<double> K2_im,
                                 const vector<vector<double> > &AllX,
                                 const vector<double> &AllWs,
                                 const vector<double> &AllNorms,
                                 double c,
                                 vector<double> Measured_X,
                                 vector<double> &ScalarProducts_re,
                                 vector<double> &ScalarProducts_im,
                                 double norm1,
                                 double norm2)
{

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Numerical integration (scalar product calculation <f2|f1> over M-cube, where M = | K1.size() - K2.size() |.
In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method - collections corresponding
to this method may be obtained using the function IntegrationCollections(...).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int N1 = K1_re.size();
    int N2 = K2_re.size();
    int M  = N2 - abs(N1 - N2); // number of variables dx_1 dx_2 ... dx_M

// ======================================================================================================================
// 1) Last step of collections of positions preparation (here we assume that N2 > N1)
// ======================================================================================================================

  /*
    vector<double> tmp;

    //X1 = AllX;

    for(unsigned i = 0; i < AllX.size(); i++)
    {
        cout << __func__ << i << " dupa " << std::endl;
        tmp = AllX.at(i);

        for (auto j : Measured_X) {
            tmp.push_back(j);
        }
//        for(unsigned j = 0; j < Measured_X.size(); j++)
//        {
//            tmp.push_back( Measured_X.at(j) );
//        }

        AddCollection(X2,tmp);
        tmp.clear();
    }

 */



// ======================================================================================================================
// 2) Summation
// ======================================================================================================================

    Complex f0, f2;
    Complex IntValue(0.,0.);

    for(unsigned s = 0; s < AllX.size(); s++)
    {
        vector<double> X1 = AllX.at(s);

        vector<double> X2, tmp;
        tmp = AllX.at(s);
        for(unsigned j = 0; j < Measured_X.size(); j++)
        {
           tmp.push_back( Measured_X.at(j) );
        }
        X2 = tmp;
        tmp.clear();

        f0 = R_WaveFunction(N1, c, X1, K1_re, K1_im, norm1);
        Complex f1( f0.real(), - f0.imag() );
        f2 = R_WaveFunction(N2, c, X2, K2_re, K2_im, norm2);
        Complex Coeffs(AllNorms.at(s) * AllWs.at(s), 0.);
        IntValue = IntValue + Coeffs * f1 * f2;

        if( s % 100000 ==0)
        {
            cout<< s << " ";
        }
    }

    ScalarProducts_re.push_back( IntValue.real() );
    ScalarProducts_im.push_back( IntValue.imag() );

    //cout<< "\n";

    //cout<< "(" << IntValue.real() << " , " << IntValue.imag() << ")" << "\n";


}













void permutator::R_ScalarProduct_approx(double P,
                                        const vector<vector<double> > &AllX,
                                        const vector<double> &AllWs,
                                        const vector<double> &AllNorms,
                                        double c,
                                        double L,
                                        vector<double> Measured_X,
                                        vector<double> &ScalarProducts_re,
                                        vector<double> &ScalarProducts_im)

{


    // ======================================================================================================================
    // 2) Summation
    // ======================================================================================================================

        Complex f0;
        Complex IntValue(0.,0.);

        Complex norm1(0.,0.);
        Complex norm2(0.,0.);

        for(unsigned s = 0; s < AllX.size(); s++)
        {
            vector<double> X1 = AllX.at(s);

            vector<double> X2, tmp;
            tmp = AllX.at(s);

            for(unsigned j = 0; j < Measured_X.size(); j++)
            {
               tmp.push_back( Measured_X.at(j) );
            }
            X2 = tmp;
            tmp.clear();


            f0 = R_pushed_WF_approx(X1, P, L, c);

            Complex f1( f0.real(), - f0.imag() );
            Complex f2( R_WF_approx(X2, L, c) , 0.);
            Complex Coeffs(AllNorms.at(s) * AllWs.at(s), 0.);

            IntValue = IntValue + Coeffs * f1 * f2;

            norm1 = norm1 + Coeffs * f1 * f0;
            norm2 = norm2 + Coeffs * f2 * f2;



            if( s % 100000 ==0)
            {
                //cout<< s << " ";
            }
        }

        ScalarProducts_re.push_back( IntValue.real() );
        ScalarProducts_im.push_back( IntValue.imag() );

        cout<< "\n\n" << "Momentum = " << P << "\n value:" << pow (fabs(IntValue), 2.) << "\n norm1: " << norm1 << "\n norm2: " << norm2;

        double sum2 = 0;

        for(unsigned y = 0; y < ScalarProducts_im.size(); y++)
        {
            double scprod = pow( ScalarProducts_im[y], 2.) + pow( ScalarProducts_re[y], 2.);

            sum2 = sum2 +scprod;
        }

        cout<< "\n\n sum of scprod^2: " << sum2;

}













void permutator::R_MC_ScalarProduct_approx(int L,
                                           double c,
                                           double P,
                                           int NumOfParticlesAfterMeasurement,
                                           int NumberOfPoints,
                                           double range,
                                           vector<double> Measured_X,
                                           vector<double> &ScalarProducts_re,
                                           vector<double> &ScalarProducts_im,
                                           vector<double> &NormsAfterMeasurement
                                           )
{// We need to choose a type of method:
//  Very easy method! Value = L^d (f(set X_1)/n + f(set X_2)/n ... f(set X_n)/n)

    double re, im, ElSize;

//========= RANGE OF INTEGRATION ======================================================================
    // I choose the range of numerical integration. In the case of very narrow maximum
    // it is enough to restrict the integration just around the maximum (soliton).
    // Of course sometimes it is... bad idea - then if the variable range > L/2 the integration
    // is over all L^d where d = NumOfParticlesAfterMeasurement.
    // Be rationall - choose range such that a > 0 and b < L.
    double a,b;
    if(range > L/2.)
    {
        a = 0;
        b = L;
        ElSize = L;
    }
    else
    {
        double NumOfMeasured = Measured_X.size();
        double SumOfMeasured = 0.;
        for(int i = 0; i < NumOfMeasured; i++)
        {
            SumOfMeasured = SumOfMeasured + Measured_X[i];
        }

        double MeanOfMeasured = SumOfMeasured/NumOfMeasured;

        a = MeanOfMeasured - range;
        b = MeanOfMeasured + range;

        ElSize = 2.* range;
    }
//====================================================================================================

//================= INTEGRATION ======================================================================
        Complex f0, f2;
        Complex Coefficient( pow( ElSize , NumOfParticlesAfterMeasurement) * 1./NumberOfPoints, 0. ); // Zmienilem z pow(L, NumOfVariables)
        Complex int_value(0.,0.);


        Complex norm1(0.,0.);
        Complex norm2(0.,0.);


        for(int s = 0; s < NumberOfPoints; s++)
        {
            //------------------------------------------------------------------------------------------------------
            vector<double> VecOfPos,VecOfPosPlusMeasured;
            for(int j = 0; j < NumOfParticlesAfterMeasurement; j++)
            {
                double RandPos =  GetRandom(a,b);
                VecOfPos.push_back( RandPos );
            }
            VecOfPosPlusMeasured = VecOfPos;
            for(int i = 0; i < Measured_X.size(); i++)
            {
                VecOfPosPlusMeasured.push_back(Measured_X[i] );
            }
            //------------------------------------------------------------------------------------------------------

            vector<double> X1 = VecOfPos;
            vector<double> X2 = VecOfPosPlusMeasured;

            f0 = R_pushed_WF_approx(X1, P, L, c);

            Complex f1( f0.real(), - f0.imag() );
            Complex f2( R_WF_approx(X2, L, c) , 0.);

            int_value = int_value +  Coefficient * f1 * f2;

            norm1 = norm1 + Coefficient * f1 * f0; // With P
            norm2 = norm2 + Coefficient * f2 * f2; // just a ground state


            if( s % 5000 ==0)
            {
                //cout<< s << " ";
            }
        }


        int_value = int_value * Complex( pow (norm1.real(),  -.5), 0. );

        re = int_value.real();
        im = int_value.imag();


    ScalarProducts_re.push_back( re );
    ScalarProducts_im.push_back( im );
    NormsAfterMeasurement.push_back(pow (norm1.real(),  -.5));


    cout<< "\n\n" << "Momentum = " << P << "\n value:" << pow (fabs(int_value), 2.) << "\n norm1: " << norm1 << "\n norm2: " << norm2;

    double sum2 = 0;

    for(unsigned y = 0; y < ScalarProducts_im.size(); y++)
    {
        double scprod = pow( ScalarProducts_im[y], 2.) + pow( ScalarProducts_re[y], 2.);

        sum2 = sum2 +scprod;
    }

    cout<< "\n\n sum of scprod^2: " << sum2;


}















// ###################################################################################################################################
// =============================== R_RelevantCoefficients ============================================================================
// ###################################################################################################################################


void permutator::R_RelevantCoefficients(vector<double> &ScalarProducts_re,
                                        vector<double> &ScalarProducts_im,
                                        vector<vector<double>> &K_Matrix_re,
                                        vector<vector<double>> &K_Matrix_im,
                                        double CutOff_Coefficient)
{
    //We need to choose only important coefficients (in fact values of scalar products calculated before)
    //CutOff_Coefficient - the relevant coefficients are not smaller than MAX[abs(coefficients)]*CutOff_Coefficient

    int n = ScalarProducts_re.size();
    int PositionOfMax = 0;
    double tmp, MaxValue;
    vector<double> AbsValues, Tmp_re, Tmp_im, Tmp_K_re, Tmp_K_im;
    vector<vector<double>> Tmp_K_Matrix_re, Tmp_K_Matrix_im;

    for(int i = 0; i < n; i++)
    {
        tmp = pow(   pow( ScalarProducts_re.at(i) , 2.) + pow( ScalarProducts_im.at(i) , 2.)   ,   0.5 );

        AbsValues.push_back( tmp );

        if( i == 0 )
        {}
        else
        {
            if( AbsValues.at( PositionOfMax ) >= tmp )
            {}
            else
            {
                PositionOfMax = i;
            }
        }
    }


    MaxValue = AbsValues.at( PositionOfMax );
    double NotLessThan = MaxValue * CutOff_Coefficient; // "real" cut off

    for(int j = 0; j < n; j++)
    {
        if( AbsValues.at(j) >= NotLessThan )
        {
            Tmp_re.push_back( ScalarProducts_re.at(j) );
            Tmp_im.push_back( ScalarProducts_im.at(j) );

            Tmp_K_re = K_Matrix_re.at(j);
            Tmp_K_im = K_Matrix_im.at(j);
            Tmp_K_Matrix_re.push_back( Tmp_K_re );
            Tmp_K_Matrix_im.push_back( Tmp_K_im );

        }
    }

    // NEW RELEVANT VALUES OF SCALAR PRODUCTS
    ScalarProducts_re = Tmp_re;
    ScalarProducts_im = Tmp_im;

    //AND CORRESPONDING K_Matrix(re and im)
    K_Matrix_re = Tmp_K_Matrix_re;
    K_Matrix_im = Tmp_K_Matrix_im;
}




// ###################################################################################################################################
// =============================== R_SCALAR_PRODUCT (CARLO) ====================================================================
// ###################################################################################################################################





void permutator::R_MonteCarlo_ScalarProduct(int Type, int L,
                                            vector<double> K1_re,
                                            vector<double> K1_im,
                                            vector<double> K2_re,
                                            vector<double> K2_im,
                                            const vector<vector<double>>& RandPos_Matrix,
                                            double c,
                                            vector<double> Measured_X,
                                            vector<double> &ScalarProducts_re,
                                            vector<double> &ScalarProducts_im,
                                            double norm1,
                                            double norm2)
{// We need to choose a type of method:
//  Type 1: Very easy method! Value = L^d (f(set X_1)/n + f(set X_2)/n ... f(set X_n)/n)
//  Type 2: Method based on the number of points counted "under the curve"

    int NumberOfPoints = RandPos_Matrix.size();
    int NumOfVariables = (RandPos_Matrix.at(0)).size();
    int N1 = K1_re.size();
    int N2 = K2_re.size();
    //here we assume that N2 > N1

    /*
    vector<vector<double>> X1, X2;
    vector<double> tmp;

    X1 = RandPos_Matrix;

    for(int i = 0; i < NumberOfPoints; i++)
    {
        tmp = RandPos_Matrix.at(i);

        for(unsigned j = 0; j < Measured_X.size(); j++)
        {
            tmp.push_back( Measured_X.at(j) );
        }

        AddCollection(X2,tmp);
        tmp.clear();
    }
    */

    double re, im;
    vector<double> X1, X2, tmp;

    if( Type == 1 )
    {
        Complex f0, f2;
        Complex Coefficient( pow(.4, NumOfVariables) * 1./NumberOfPoints, 0. ); // Zmienilem z pow(L, NumOfVariables)
        Complex int_value(0.,0.);

        for(int s = 0; s < NumberOfPoints; s++)
        {
            X1 = RandPos_Matrix.at(s);
            tmp = RandPos_Matrix.at(s);

            for(unsigned j = 0; j < Measured_X.size(); j++)
            {
                tmp.push_back( Measured_X.at(j) );
            }
            X2 = tmp;

            f0 = R_WaveFunction(N1, c, X1, K1_re, K1_im, norm1);
            Complex f1( f0.real(), - f0.imag() );
            f2 = R_WaveFunction(N2, c, X2, K2_re, K2_im, norm2);

            int_value = int_value +  Coefficient * f1 * f2;

            if( s % 5000 ==0)
            {
                //cout<< s << " ";
            }
        }

        re = int_value.real();
        im = int_value.imag();


    }
    else{

    Complex f0, f2, value;
    double ValRe, ValIm, RandRe, RandIm;
    double range = 3.;
    int NumRe = 0;
    int NumIm = 0;

        for(int s = 0; s < X1.size(); s++)
        {
            X1 = RandPos_Matrix.at(s);
            tmp = RandPos_Matrix.at(s);

            for(unsigned j = 0; j < Measured_X.size(); j++)
            {
                tmp.push_back( Measured_X.at(j) );
            }
            X2 = tmp;

            f0 = R_WaveFunction(N1, c, X1, K1_re, K1_im, norm1);
            Complex f1( f0.real(), - f0.imag() );
            f2 = R_WaveFunction(N2, c, X2, K2_re, K2_im, norm2);

            value = f1 * f2;
            ValRe = value.real();
            ValIm = value.imag();

            RandRe = GetRandom( -range , range);
            RandIm = GetRandom( -range , range);

//---------------- Re part -------------------------------------------------
            if( ValRe > 0. )
            {
                if( RandRe > 0. )
                {
                    if( ValRe > RandRe )
                    {
                        NumRe = NumRe + 1;
                    }
                }
            }
            else
            {
                if( RandRe < 0. )
                {
                    if( ValRe < RandRe )
                    {
                        NumRe = NumRe - 1;
                    }
                }
            }

//---------------- Im part -------------------------------------------------
            if( ValIm > 0. )
            {
                if( RandIm > 0. )
                {
                    if( ValIm > RandIm )
                    {
                        NumIm = NumIm + 1;
                    }
                }
            }
            else
            {
                if( RandIm < 0. )
                {
                    if( ValIm < RandIm )
                    {
                        NumIm = NumIm - 1;
                    }
                }
            }


        }

        re = 2. * range * pow( L ,NumOfVariables ) * (1.* NumRe)/NumberOfPoints;
        im = 2. * range * pow( L ,NumOfVariables ) * (1.* NumIm)/NumberOfPoints;

    }



    ScalarProducts_re.push_back( re );
    ScalarProducts_im.push_back( im );


}







// ###################################################################################################################################
// ============================ R_WaveFunctionAfterMeasurement =======================================================================
// ###################################################################################################################################


std::complex<double> permutator::R_WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix_re,
                                                                vector<vector<double> > K_Matrix_im,
                                                                vector<double> ScalarProducts_re,
                                                                vector<double> ScalarProducts_im,
                                                                vector<double> Remaining_X,
                                                                double c,
                                                                double phi,
                                                                double t,
                                                                vector<double> Norms)
{
// The wave function after measurement of particle positions - it is given by a sum over N-M particle WaveFunctions with proper coefficients (scalar products)
// Cause of the structure of the Wave Function we obtain time evolution in very easy way.

    int n = Remaining_X.size();
    int NumOfRows = K_Matrix_re.size();

    double re = 0., im = 0.;
    int i;

//#pragma omp parallel  private(i) reduction(+:re) reduction(+:im)
//{
    Complex WF(0. , 0.);
//    #pragma omp for

    for(i = 0; i < NumOfRows; i++)
    {
        double energy = 0.;
        double momentum = 0.;

        for(int j = 0; j < n; j++)
        {
            energy = energy + pow( (K_Matrix_re.at(i)).at(j) , 2.) - pow( (K_Matrix_im.at(i)).at(j) , 2.); // Cause terms 2 I Re(K^2) Im(K^2) should cancel out
            momentum = momentum + (K_Matrix_re.at(i)).at(j);
        }

        Complex Coefficient( ScalarProducts_re.at(i) , ScalarProducts_im.at(i) );
        Complex I_Exp_t(0. , (energy + phi * momentum) * t);
        Complex EvolutionExp;
        EvolutionExp = exp(- I_Exp_t);

        WF = WF +  EvolutionExp * Coefficient  * R_WaveFunction(n, c, Remaining_X, K_Matrix_re.at( i ), K_Matrix_im.at( i ), Norms.at(i) );

        //cout<< "\n time t=" << t << "\n\n";
        //cout<< "\n evolution exp =" << "(" << EvolutionExp.real() << " , " <<  EvolutionExp.imag() <<  ")" << "\n\n";

    }

    re=WF.real();
    im=WF.imag();
//}

Complex Wavefunction(re,im);

    return Wavefunction;

}




// ###################################################################################################################################
// ============================== R_Metropolis_TimeDependent =======================================================================
// ###################################################################################################################################


void permutator::R_Metropolis_TimeDependent(int N,
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
                                            vector<double> Norms)
{
    double ratio;
    double dif;
    double random0to1;
    double wf0;
    double wf1;

    for(int i = 0; i < N; i++)
    {
        Positions(X0, X1, delta, i, L);
    }
        wf0 =  pow( fabs( R_WaveFunctionAfterMeasurement(K_Matrix_re, K_Matrix_im, ScalarProducts_re, ScalarProducts_im, X0, c, phi, t, Norms)  ), 2.);
        wf1 =  pow( fabs( R_WaveFunctionAfterMeasurement(K_Matrix_re, K_Matrix_im, ScalarProducts_re, ScalarProducts_im, X1, c, phi, t, Norms)  ), 2.);

        ratio = wf1/wf0;

        if( ratio >= 1 )
        {
            WaveFValues.push_back(wf1);
            AddCollection(Chain,X1);
        }
        else
        {
            random0to1 = GetRandom(0,1);
            dif = ratio - random0to1;
                if( dif > 0 )
                {
                    WaveFValues.push_back(wf1);
                    AddCollection(Chain,X0);
                }
                else
                {
                    X1 = X0;
                }
        }
     X0 = X1;
}


permutator _Permutator;













