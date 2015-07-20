#ifndef SPERMUTATOR_H
#define SPERMUTATOR_H

#include <QObject>

class sPermutator : public QObject
{
    Q_OBJECT
public:
    explicit sPermutator(QObject *parent = 0);
//    ~sPermutator();


private:
    QList<QList<int> > sperms;
    int _N = 0;
    QList<int> _sPermPool;

public:
    void setSpermSize(int N);
    int  getSpermSize();
    void sPermSwap(int a, int b);
    void sPermute(int k, int size);
    void appendSPerm(int size);

signals:

public slots:
};

#endif // SPERMUTATOR_H
