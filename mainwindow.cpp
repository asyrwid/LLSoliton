#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
// initial values!
//*************************************************************************
//*************************************************************************

    int    N = 4; // number of particles
    double c = 1; // coupling constant
    double L = 1; // system length
    double delta = 0.1; // maximal "jump" between single particle
                        // position from one realization to another
    QList<double> X0 = {.1,.3,.5,.7}; // initial collection of particle positions
    QList<double> K = {-2.24633, -0.714948, 0.714948, 2.24633}; //quasimomentas
    QList<double> X1 = X0; //
//*************************************************************************
//*************************************************************************



/*
//N! calculation which we need to know how many permutations we will have
//*************************************************************************
        int Nfact = 1; //Nfact = N! (Nfactorial)
        for (int i = 2; i <= N; i++) //Nfactorial value (N!) calculation
        { Nfact = Nfact * i; }
        //qDebug() << Nfact;
// Actually, it is the number of rows of "perms"
// *************************************************************************
*/



    ui->setupUi(this);
    //qDebug() << _sPermutator.getSpermSize(); // initially _N=0;
    _sPermutator.setSpermSize(N); // de facto particle number
    //qDebug() << _sPermutator.getSpermSize(); // N is established now!
    _sPermutator.sPermute(N,N); // all permutations preparation
    //_sPermutator.sProduct(N,0, c, X,K);
    //qDebug() << _sPermutator.SignFunction(1.2, 2.3);
    //_sPermutator.WaveFunction(N,c,X,K);
    _sPermutator.Positions(X0,X1,delta,1,L);
    _sPermutator.Positions(X0,X1,delta,1,L);

}

MainWindow::~MainWindow()
{
    delete ui;
}
