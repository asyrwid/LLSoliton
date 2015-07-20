#include "spermutator.h"
#include <QDebug>
sPermutator::sPermutator(QObject *parent) : QObject(parent)
{

}

void sPermutator::setSpermSize(int N)
{
    _N = N;
    _sPermPool.clear();
    for (int i=0; i<N; i++)
    {
        _sPermPool.append(i+1);
    }
}

int sPermutator::getSpermSize()
{
    return _N;
}

void sPermutator::sPermSwap(int a, int b)
{
    if (a == b)
    {
     return;
    }

    _sPermPool.swap(a,b);

    _sPermity = !_sPermity;
}

void sPermutator::appendSPerm(int size)
{
    sperms.append(_sPermPool.mid(0, size));
    if (_sPermity == true){
        _sPermSigns.append(1);
    }
    else{
        _sPermSigns.append(-1);
    }
    qDebug() << sperms.last() << ' ' << _sPermSigns.last();
}

void sPermutator::sPermute(int k, int size)
{
    //qDebug() << k << __func__;
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

