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
{                               //the particle number
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


//*******************************************************************************************************************
//*******************************************************************************************************************
//WAVE FUNCTION PREPARATION

typedef std::complex<double> Complex;

std::complex<double> permutator::sProduct(int n, int PermNumber, double c, vector<double> X, vector<double> K)
{   // 1.
    //Product (j>k which appears in the wave function) definition
    //PermNumber is the number of permutation in permutations array (row)
    Complex  prod(1.,0.);
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < j; k++)
        {
        prod = prod* Complex( K[GetPermEl(PermNumber,j)] - K[GetPermEl(PermNumber,k)], -c*SignFunction( X[j], X[k]) ) ;
        }

    }
    // cout << prod.real() << ',' << prod.imag();

    //2.
    //Exponent with a sum over X[]K[]
    double sumexp = 0.;
    for(int j = 0; j < n; j++)
    {
        sumexp = sumexp + K[GetPermEl(PermNumber,j)] * X[j];
    }
    Complex ImSumExp(0., sumexp);
    Complex ImExp;
    ImExp = exp(ImSumExp);
    // cout << ImExp.real() << ',' << ImExp.imag();

    //3.
    //sign of the permutation given by PermNumber
    int SPV = _PermSigns.at(PermNumber); // + or - 1 (Parity)
    Complex SignValue(SPV, 0);

    //4.
    //Finally, we can calculate 1 element of "big" sum apearing in the wave function
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
    //N! calculation which we need to know how many permutations we will have
    //*************************************************************************
            int Nfact = 1; //Nfact = N! (Nfactorial)
            for (int i = 2; i <= n; i++) //Nfactorial value (N!) calculation
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
    Complex wavefunction(re,im);

    //Normalization factor

    double norm = Nfact;
    for (int j = 0; j < n ; j++)
    {
        for (int k = 0; k < j; k++)
        {
        norm = norm * ( pow ((K[j]-K[k]),2.) + pow (c,2.) );
        }
    }

    Complex Norm( 1./(sqrt( norm)),0 ) ;
    wavefunction = Norm * wavefunction;
   // cout << abs(wavefunction) <<  "\n";

    //cout << wavefunction.real() << ',' << wavefunction.imag();

    //wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return wavefunction;
}

//*******************************************************************************************************************
//*******************************************************************************************************************

//COLLECTIONS OF PARTICLE POSITIONS GENERATION
//*******************************************************************************************************************
//*******************************************************************************************************************

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
        wf0 =  pow (abs(WaveFunction(N,c, X0, K) ), 2.);
        wf1 =  pow (abs(WaveFunction(N,c, X1, K) ), 2.);

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


// Finding the minimum of "the last" probability density
// we need to find a jump(s) of the phase of wave function
//*******************************************************************************************************************
//*******************************************************************************************************************

vector<double> permutator::PhaseCalculation(std::complex<double> WF)
{// gives a vector (abs(Wavefunction), Phase), value of Wave Function is needed to find appropriate jump

    double PI_Value = 3.14159265358979323846;
    double  WF_re, WF_im, Phase;
    vector<double> Out;
    WF_re = WF.real();
    WF_im = WF.imag();

    if( WF_re > 0){
        if( WF_im > 0 )
        {
            Phase = atan( WF_im/WF_re );
        }
        else{// means WF_re > 0, WF_im <0
            Phase = PI_Value - atan( abs( WF_im/WF_re ) );
        }
    }
    else{// means WF_re < 0
        if( WF_im > 0 )
        {
            Phase = 2.*PI_Value - atan( abs( WF_im/WF_re ) );
        }
        else{// means WF_re < 0, WF_im < 0
            Phase = PI_Value + atan( abs( WF_im/WF_re ) );
        }
    }

    Out.push_back( abs(WF) );
    Out.push_back( Phase );

    return Out;
}


void permutator::PhaseJump(vector<double> &Jumps, int InitialDivision, int n, int ChainElement, double JumpRestriction ,double c, double L, vector<vector<double> > Chain, vector<double> K)
{
    vector<double>  PosX = Chain.at(ChainElement); //Positions from the ChainElement of Chain
    vector<double>  InitialPositions, PhaseValues, WFValues, temp;
    double  Phase, WF;

//Preparation of initial set of points (0,L*1/(InitialDiv-1),...,L)
    for(int i = 0; i < InitialDivision; i++)
    {
        InitialPositions.push_back( i*L/(InitialDivision - 1) );
    }

//Phase and abs of Wave Function in all the points from initial set
    for(int i = 0; i < InitialDivision; i++)
    {
        PosX.at(0) = InitialPositions.at(i);
        temp = PhaseCalculation(WaveFunction(n, c, PosX, K) );
        Phase = temp.at(1);
        WF = temp.at(0);
        PhaseValues.push_back( Phase );
        WFValues.push_back( WF );

        cout << i << "   " << PhaseValues.at(i) << "   " << WFValues.at(i) <<"\n";
    }

//Jump searching
    vector<double> JumpValue;
    vector<int> JumpPositions, JumpPositions2;
    double diff;

    for(int i = 0; i < InitialDivision-1; i++)
    {
        diff = sqrt( pow ( PhaseValues.at(i+1) - PhaseValues.at(i), 2) );
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

    int s, p;
    for(int j = 0; j < JumpPositions.size(); j++)
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
    cout << "\n" << "============" << "\n";
    cout << p << "\n";

    Jumps = JumpValue;


}















