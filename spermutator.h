#ifndef SPERMUTATOR_H
#define SPERMUTATOR_H
#include <complex>
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
    QList<int> _sPermSigns;
    bool _sPermity = true;

public:
    void setSpermSize(int N);
    int  getSpermSize();
    void sPermSwap(int a, int b);
    void sPermute(int k, int size);
    void appendSPerm(int size);
    std::complex<double> sProduct(int n, int PermNumber, double c, QList<double> X, QList<double> K);
    int SignFunction(double x1, double x2);
    int GetPermEl(int i, int j);
    std::complex<double> WaveFunction(int n, double c, QList<double> X, QList<double> K);
    double GetRandom(double min, double max);
    QList<double> ParticlePositions;
    void Positions(QList<double> X0, QList<double> &X1, double delta, int n, double L);
    void AddCollection(QList<QList<double> > &Chain, QList<double>  Set);
    void Metropolis(int N,double c, QList<QList<double> > &Chain,  QList<double> &WaveFValues, QList<double>  &X0, QList<double>  &X1, QList<double>  K, double delta, int n, double L);

signals:

public slots:
};

#endif // SPERMUTATOR_H
