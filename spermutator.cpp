#include "spermutator.h"
#include <QDebug>
#include <iterator>
#include <math.h>
#include <random>

sPermutator::sPermutator(QObject *parent) : QObject(parent)
{
}

void sPermutator::setSpermSize(int N) //particle number setting
{                                     //and initial set 1,...,N
    _N = N;                           //preparation
    _sPermPool.clear();
    for (int i=0; i<N; i++)
    {
        _sPermPool.append(i);
    }
}

int sPermutator::getSpermSize() //expectorate (wykrztuś,wypluj)
{                               //the particle number
    return _N;
}

void sPermutator::sPermSwap(int a, int b) //"swap" additional operator
{
    if (a == b)
    {
     return;
    }

    _sPermPool.swap(a,b);

    _sPermity = !_sPermity; //de facto Parity change after
}                           //swap operation :)


void sPermutator::appendSPerm(int size) //operation appending
{                                //another permutation + sign
    sperms.append(_sPermPool.mid(0, size));
    //qDebug() << _sPermPool;
    if (_sPermity == true){
        _sPermSigns.append(1);
    }
    else{
        _sPermSigns.append(-1);
    }
    //qDebug() << _sPermSigns.last();
    //qDebug() << sperms.last() << ' ' << _sPermSigns.last();
}

void sPermutator::sPermute(int k, int size)//all permutations
{                        //preparation + signs (appendSPerm used)
    if(k == 0)
    {
        appendSPerm(size);
    }
    else{
       for(int i = k-1; i >= 0; i--)
            {
            sPermSwap(i,k-1);
            sPermute(k-1,size);
            sPermSwap(i,k-1);
            }
    }
}

int sPermutator::GetPermEl(int i, int j) // acces to (i,j) element of
{                                        // all permutations array
    return (sperms.at(i)).at(j);
}


int sPermutator::SignFunction(double x1, double x2) // sign function
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

std::complex<double> sPermutator::sProduct(int n, int PermNumber, double c, QList<double> X, QList<double> K)
{   // 1.
    //Product (j>k which appears in the wave function) definition
    //PermNumber is the number of permutation in permutations array (row)
    Complex  prod(1.,0.);
    for (int j = 0; j < n  ; j++)
    {
        for (int k = 0; k < j; k++)
        {
        prod = prod* Complex( K[GetPermEl(PermNumber,j)] - K[GetPermEl(PermNumber,k)], -c*SignFunction( X[j], X[k]) ) ;
        }

    }
    // qDebug() << prod.real() << ',' << prod.imag();

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
    // qDebug() << ImExp.real() << ',' << ImExp.imag();

    //3.
    //sign of the permutation given by PermNumber
    int SPV = _sPermSigns.at(PermNumber); // + or - 1 (Parity)
    Complex SignValue(SPV, 0);
    //4.
    //Finally, we can calculate 1 element of "big" sum apearing
    //in the wave function
    Complex ProductElement;
    ProductElement = SignValue * ImExp * prod;
/*
    qDebug() <<'('<<ProductElement.real()<<','<<ProductElement.imag()<<')' ;
    qDebug() <<'('<<SignValue.real()<<','<<SignValue.imag()<<')' ;
    qDebug() << prod.real() << ImExp.real();
    qDebug() << prod.imag() << ImExp.imag();
*/

   // qDebug() << ProductElement.real() << ',' << ProductElement.imag();
    return ProductElement;
}

std::complex<double> sPermutator::WaveFunction(int n, double c, QList<double> X, QList<double> K)
{
    //N! calculation which we need to know how many permutations we will have
    //*************************************************************************
            int Nfact = 1; //Nfact = N! (Nfactorial)
            for (int i = 2; i <= n; i++) //Nfactorial value (N!) calculation
            { Nfact = Nfact * i; }
            //qDebug() << Nfact;
    // Actually, it is the number of rows of "perms"
    // ********r****************************************************************
    double re=0,im=0;
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
        norm = norm * ( (K[j]-K[k])*(K[j]-K[k]) + c*c );
        }
    }
    Complex Norm( 1/(sqrt(norm)),0 ) ;
    wavefunction = Norm * wavefunction;

    //qDebug() << wavefunction.real() << ',' << wavefunction.imag();

    wavefunction = (wavefunction.real())*(wavefunction.real()) + ( wavefunction.imag())*( wavefunction.imag()) ;
    return wavefunction;
}

//*******************************************************************************************************************
//*******************************************************************************************************************

//COLLECTIONS OF PARTICLE POSITIONS GENERATION
//*******************************************************************************************************************
//*******************************************************************************************************************

double sPermutator::GetRandom(double min, double max)
{// Returns a random double between min and max
    return ((double) rand()*(max-min)/(double)RAND_MAX + min);
}
void sPermutator::Positions(QList<double> X0, QList<double> &X1, double delta, int n, double L)
{// one step in Metropolis x1(n)->x1'(n) = x0(n)+ delta*Q, Qin [-L,L]
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

    //qDebug() << X1;
}


void sPermutator::AddCollection(QList<QList<double> > &Chain, QList<double>  Set)
{
    Chain.append( Set.mid(0, Set.size() ) );
    //qDebug() << Chain;
}

void sPermutator::Metropolis(int N,double c, QList<QList<double> > &Chain, QList<double> &WaveFValues,  QList<double>  &X0,QList<double>  &X1, QList<double>  K, double delta, int n, double L)
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
        wf0 =  WaveFunction(N,c, X0, K).real();
        wf1 =  WaveFunction(N,c, X1, K).real();

        ratio = wf1/wf0;

        if( ratio >= 1 )
        {
            WaveFValues.append(wf1);
            AddCollection(Chain,X1);
        }
        else
        {
            random0to1 = GetRandom(0,1);
            dif = ratio - random0to1;
                if( dif > 0 )
                {
                    WaveFValues.append(wf1);
                    AddCollection(Chain,X0);
                }
                else
                {
                    X1 = X0;
                }
        }
     X0 = X1;

}











