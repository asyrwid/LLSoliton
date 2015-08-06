#include "permutator.h"
#include <iterator>
#include <math.h>
#include <random>

//permutator::permutator()
//{
//`}

void permutator::setPermSize(int N) //particle number setting
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
}                           //swap operation :)


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
    if(x1 < x2){
        t = -1;}
    else{
        t = 1;}
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
    Complex SignValue(SPV, 0);

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

std::complex<double> permutator::WaveFunction(int n, double c, vector<double> X, vector<double> K)
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
            wavefunction = wavefunction + sProduct(n, i, c, X, K);
        }
        re=wavefunction.real();
        im=wavefunction.imag();
    }
    Complex Wavefunction(re,im);

    // Normalization factor

    double norm = Nfact;
    for (int j = 0; j < n ; j++)
    {
        for (int k = 0; k < j; k++)
        {
        norm = norm * ( pow ((K[j]-K[k]),2.) + pow (c,2.) );
        }
    }

    Complex Norm( 1./(sqrt( norm)),0 ) ;
    Wavefunction = Norm * Wavefunction;
    //cout << fabs(wavefunction) <<  "\n";

    //cout << wavefunction.real() << ',' << wavefunction.imag();

    //wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return Wavefunction;
}

//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------





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

void permutator::Metropolis(int N,double c, vector<vector<double> > &Chain, vector<double> &WaveFValues,  vector<double>  &X0,vector<double>  &X1, vector<double>  K, double delta, int n, double L)
{
    double ratio;
    double dif;
    double random0to1;
    double wf0;
    double wf1;

    for(int i = 0; i < n; i++)
    {
        Positions(X0, X1, delta, i, L);
    }
        wf0 =  pow (fabs(WaveFunction(N,c, X0, K) ), 2.);
        wf1 =  pow (fabs(WaveFunction(N,c, X1, K) ), 2.);

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





//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------
//***************************** FINDING THE MINIMUM OF "THE LAST" PROBABILITY DENSITY ************************************************
//**************************** WE NEED TO FIND A JUMP(S) OF THE PHASE OF WAVE FUNCTION ***********************************************
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

vector<double> permutator::PhaseCalculation(std::complex<double> WF)
{// gives a vector (fabs(Wavefunction), Phase), value of Wave Function is needed to find appropriate jump
 // we can use collected values to plot Prob density and phase distribution

    double  Phase;
    vector<double> Out;

    Phase = arg ( WF );

    Out.push_back( pow ( pow(10., - 150.) * fabs (WF), 2 ) ); // Prob density value
    Out.push_back( Phase ); // Phase value

    return Out;
}




void permutator::MinimaFinder(vector<double> &Minima, vector<vector<double>> &PhaseMatrix, vector<vector<double>> &WFMatrix, int NumOfBisections, int InitialDivision, int n, int ChainElement, double JumpRestriction ,double c, double L, vector<vector<double> > Chain, vector<double> K)
{
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The MininaFinder function find the value of proper minimum of prob density based on the phase jump and value of prob
density in the point corresponding to the phase discontinuity. We know that the jump phase shouldn't be bigger than Pi :)

    Finding the values of phase and prob density in the InitialDivision number of points we can find an interesting
bracket. After that we use the bisection method. Phase and Prob density plot points obtained during initial calculation
may be used to prepare plots - so we can kill two birds with one stone.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


// ======================================================================================================================
// 1) Calculation of the Phase and Prob density values in initial points
// ======================================================================================================================

    vector<double>  PosX = Chain.at(ChainElement); // Choosing the collection of positions from the Chain (  .at(ChainElement)  )
    vector<double>  InitialPositions, PhaseValues, WFValues, temp;
    double  Phase, WF;

//Preparation of initial set of points (0,L*1/(InitialDiv-1),...,L)

    for(int i = 0; i < InitialDivision; i++)
    {
        InitialPositions.push_back( i*L/(InitialDivision - 1) );
    }

//Values of Phase and fabs(Wave Function)^2 in all the points from initial set

    for(int i = 0; i < InitialDivision; i++)
    {
        PosX.at(0) = InitialPositions.at(i);
        temp = PhaseCalculation(WaveFunction(n, c, PosX, K) );

        Phase = temp.at(1);
        WF = temp.at(0);

        PhaseValues.push_back( Phase ); // Set of Phase values in initial points
        WFValues.push_back( WF );       // Set of Prob denisty values in initial points

       // cout << i << "   " << PhaseValues.at(i) << "   " << WFValues.at(i) <<"\n";
    }

//Now we want to append values of Phase and Prob density to matrices


    AddCollection(PhaseMatrix, PhaseValues);
    AddCollection(WFMatrix, WFValues);


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
    WFd = fabs ( WaveFunction(n, c, PosX, K) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction(n, c, PosX, K) );



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
    WFd = fabs ( WaveFunction(n, c, PosX, K) );
    PosX.at(0) = E;
    WFe = fabs ( WaveFunction(n, c, PosX, K) );
}
  //  cout << "\n" << B << "   " << WFb << "\n";
    Minima.push_back( B );

}












