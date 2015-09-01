#include "permutator.h"
#include <iterator>
#include <math.h>
#include <random>
#include <omp.h>

//permutator::permutator()
//{
//`}



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



std::complex<double> permutator::sProduct(int n, int PermNumber, double c, vector<double> X, vector<double> K)
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








std::complex<double> permutator::WaveFunction(int n, double c, vector<double> X, vector<double> K, double norm)
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


    Complex Norm( pow(norm2, - 0.5) * pow(norm, - 0.5) ,0. ) ;
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







void permutator::Positions(vector<double> X0, vector<double> &X1, double delta, int n, double L)
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










void permutator::AddCollection(vector<vector<double> > &Chain, vector<double>  Set)
{
    Chain.push_back( Set );
    //cout << Chain;
}








void permutator::Metropolis(int N,double c, vector<vector<double> > &Chain, vector<double> &WaveFValues,  vector<double>  &X0,vector<double>  &X1, vector<double>  K, double delta, double L, double norm)
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
        wf0 =  pow (fabs(WaveFunction(N,c, X0, K, norm) ), 2.);
        wf1 =  pow (fabs(WaveFunction(N,c, X1, K, norm) ), 2.);

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

    Out.push_back( pow ( pow(10., - 150.) * fabs (WF), 2 ) ); // Prob density value (position 0)
    Out.push_back( Phase ); // Phase value (position 1)
    Out.push_back( Sin ); // Sin value (position 2)
    Out.push_back( Cos ); // Cos value (position 3)

    return Out;
}





void permutator::MinimaFinder(vector<double> &Minima, vector<vector<double>> &PhaseMatrix, vector<vector<double>> &WFMatrix, vector<vector<double> > &SinMatrix, vector<vector<double> > &CosMatrix, int NumOfBisections, int InitialDivision, int n, int ChainElement, double JumpRestriction ,double c, double L, vector<vector<double> > Chain, vector<double> K, double norm)
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





//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//***************************** SCALAR PRODUCT - NUMERICAL INTEGRATION ***************************************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------


void permutator::IntegrationCollections(vector<double> K1, vector<double> K2, vector<vector<double>> &AllX, vector<double> &AllWs, vector<double> &AllNorms, double L, int PolynomialDegree, int IntervalDivision)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Preparation of collections which are indispensable to numerical integration (scalar product calculation <f2|f1>
over M-cube, where M = | K1.size() - K2.size() |.

Of course Measured_X.size() + M = N

In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int N1 = K1.size();
    int N2 = K2.size();
    int M  = abs(N1 - N2); // number of variables dx_1 dx_2 ... dx_M

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
                CombSize1 = pow ( tmpAllX.size(), - 2. + m );
                CombSize2 = pow ( tmpAllX.size(), - 1. + m ) + 1;
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











void permutator::ScalarProduct(vector<double> K1, vector<double> K2, vector<vector<double>> AllX, vector<double> AllWs, vector<double> AllNorms, double c, vector<double> Measured_X, vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im, double norm1, double norm2)
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


    //#pragma omp for
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
    //double c2 = (c * N1 )/N2;

//#pragma omp parallel  private(s) reduction(+:re) reduction(+:im)
//{
    Complex int_value(0.,0.);

//    #pragma omp for
    for(s = 0; s < X1.size(); s++)
    {
        f0 = WaveFunction(N1, c, X1.at(s), K1, norm1);
        Complex f1( f0.real(), - f0.imag() );
        f2 = WaveFunction(N2, c, X2.at(s), K2, norm2);

        int_value = int_value + AllNorms.at(s) * AllWs.at(s) * f1 * f2;

        if( s % 500 ==0)
        {
            cout<< s << " ";
        }
    }

    re=int_value.real();
    im=int_value.imag();

//}

    Complex IntValue(re,im);

    ScalarProducts_re.push_back( IntValue.real() );
    ScalarProducts_im.push_back( IntValue.imag() );
    //cout<< "\n";

    //cout<< "(" << IntValue.real() << " , " << IntValue.imag() << ")" << "\n";


}





//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//***************** RELEVANT COEFFICIENTS AND WAVE FUNCTIONS AFTER MEASUREMENT *******************************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

void permutator::RelevantCoefficients(vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im, vector<vector<double>> &K_Matrix, double CutOff_Coefficient, vector<double> &Norms)
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







std::complex<double> permutator::WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, vector<double> Remaining_X, double c, double t, vector<double> Norms)
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







void permutator::Metropolis_TimeDependent(int N,double c, double t, vector<vector<double> > &Chain, vector<double> &WaveFValues,  vector<double>  &X0,vector<double>  &X1, vector<vector<double>>  K_Matrix, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, double delta, double L, vector<double> Norms)
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
        wf0 =  pow( fabs( WaveFunctionAfterMeasurement(K_Matrix, ScalarProducts_re, ScalarProducts_im, X0, c, t, Norms)  ), 2.);
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
                }
        }
     X0 = X1;
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


std::complex<double> permutator::R_sProduct(int n, int PermNumber, double c, vector<double> X, vector<double> K_re, vector<double> K_im)
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

std::complex<double> permutator::R_WaveFunction(int n, double c, vector<double> X, vector<double> K_re, vector<double> K_im)
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
    #pragma omp parallel  private(i) reduction(+:re) reduction(+:im)
    {
        Complex wavefunction(0.,0.);
        #pragma omp for
        for(i = 0; i < Nfact; i++)
        {
            wavefunction = wavefunction + R_sProduct(n, i, c, X, K_re, K_im);
        }
        re = wavefunction.real();
        im = wavefunction.imag();
    }
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

    Wavefunction =  Norm *  Wavefunction;
    //cout << fabs(wavefunction) <<  "\n";
    //cout << Norm.real() << ", " << Norm.imag() << "\n";
    //cout << Wavefunction.real() << ',' << Wavefunction.imag() << "\n";

    //wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return Wavefunction;
}



// ###################################################################################################################################
// =============================== R_IntegrationCollections ==========================================================================
// ###################################################################################################################################


void permutator::R_IntegrationCollections(vector<double> K1_re, vector<double> K2_re, vector<vector<double>> &AllX, vector<double> &AllWs, vector<double> &AllNorms, double L, int PolynomialDegree, int IntervalDivision)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Preparation of collections which are indispensable to numerical integration (scalar product calculation <f2|f1>
over M-cube, where M = | K1.size() - K2.size() |.

Of course Measured_X.size() + M = N

In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int N1 = K1_re.size();
    int N2 = K2_re.size();
    int M  = abs(N1 - N2); // number of variables dx_1 dx_2 ... dx_M

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
                CombSize1 = pow ( tmpAllX.size(), - 2. + m );
                CombSize2 = pow ( tmpAllX.size(), - 1. + m ) + 1;
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


void permutator::R_ScalarProduct(vector<double> K1_re, vector<double> K1_im, vector<double> K2_re, vector<double> K2_im, vector<vector<double>> AllX, vector<double> AllWs, vector<double> AllNorms, double c, vector<double> Measured_X, vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Numerical integration (scalar product calculation <f2|f1> over M-cube, where M = | K1.size() - K2.size() |.
In order to perform the calculation we use Gauss-Legendre quadrature multidimensional method - collections corresponding
to this method may be obtained using the function IntegrationCollections(...).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    int N1 = K1_re.size();
    int N2 = K2_re.size();

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
    Complex IntValue(0.,0.);

    for(unsigned s = 0; s < X1.size(); s++)
    {
        f0 = R_WaveFunction(N1, c, X1.at(s), K1_re, K1_im);
        Complex f1( f0.real(), - f0.imag() );
        f2 = R_WaveFunction(N2, c, X2.at(s), K2_re, K2_im);

        IntValue = IntValue + AllNorms.at(s) * AllWs.at(s) * f1 * f2;

        if( s % 500 ==0)
        {
            cout<< s << " ";
        }
    }

    ScalarProducts_re.push_back( IntValue.real() );
    ScalarProducts_im.push_back( IntValue.imag() );

    //cout<< "\n";

    //cout<< "(" << IntValue.real() << " , " << IntValue.imag() << ")" << "\n";


}



// ###################################################################################################################################
// =============================== R_RelevantCoefficients ============================================================================
// ###################################################################################################################################


void permutator::R_RelevantCoefficients(vector<double> &ScalarProducts_re, vector<double> &ScalarProducts_im, vector<vector<double>> &K_Matrix_re, vector<vector<double>> &K_Matrix_im, double CutOff_Coefficient)
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
// ============================ R_WaveFunctionAfterMeasurement =======================================================================
// ###################################################################################################################################


std::complex<double> permutator::R_WaveFunctionAfterMeasurement(vector<vector<double> > K_Matrix_re, vector<vector<double> > K_Matrix_im, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, vector<double> Remaining_X, double c, double t)
{
// The wave function after measurement of particle positions - it is given by a sum over N-M particle WaveFunctions with proper coefficients (scalar products)
// Cause of the structure of the Wave Function we obtain time evolution in very easy way.

    int n = Remaining_X.size();
    int NumOfRows = K_Matrix_re.size();

    double re = 0., im = 0.;
    int i;

#pragma omp parallel  private(i) reduction(+:re) reduction(+:im)
{
    Complex WF(0. , 0.);
    #pragma omp for

    for(i = 0; i < NumOfRows; i++)
    {
        double energy = 0;

        for(int j = 0; j < n; j++)
        {
            energy = energy + pow( (K_Matrix_re.at(i)).at(j) , 2.) - pow( (K_Matrix_im.at(i)).at(j) , 2.);
        }

        Complex Coefficient( ScalarProducts_re.at(i) , ScalarProducts_im.at(i) );
        Complex I_Exp_t(0. , energy * t);
        Complex EvolutionExp;
        EvolutionExp = exp(- I_Exp_t);

        WF = WF +  EvolutionExp * Coefficient  * R_WaveFunction(n, c, Remaining_X, K_Matrix_re.at( i ), K_Matrix_im.at( i ) );

        //cout<< "\n time t=" << t << "\n\n";
        //cout<< "\n evolution exp =" << "(" << EvolutionExp.real() << " , " <<  EvolutionExp.imag() <<  ")" << "\n\n";

    }

    re=WF.real();
    im=WF.imag();
}

Complex Wavefunction(re,im);

    return Wavefunction;

}




// ###################################################################################################################################
// ============================== R_Metropolis_TimeDependent =======================================================================
// ###################################################################################################################################


void permutator::R_Metropolis_TimeDependent(int N,double c, double t, vector<vector<double> > &Chain, vector<double> &WaveFValues,  vector<double>  &X0,vector<double>  &X1, vector<vector<double>>  K_Matrix_re, vector<vector<double>>  K_Matrix_im, vector<double> ScalarProducts_re, vector<double> ScalarProducts_im, double delta, double L)
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
        wf0 =  pow( fabs( R_WaveFunctionAfterMeasurement(K_Matrix_re, K_Matrix_im, ScalarProducts_re, ScalarProducts_im, X0, c, t)  ), 2.);
        wf1 =  pow( fabs( R_WaveFunctionAfterMeasurement(K_Matrix_re, K_Matrix_im, ScalarProducts_re, ScalarProducts_im, X1, c, t)  ), 2.);

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








