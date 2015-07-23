#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QCloseEvent>
#include <QTableWidget>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
// initial values!
//*************************************************************************
//*************************************************************************

    int    N = 9; // number of particles
    double c = 1; // coupling constant
    double L = 1; // system length
    double delta = 0.2; // maximal "jump" between single particle
                        // position from one realization to another

    // initial collection of particle positions
    QList<double> X0 = {.1,.3,.5,.6,.1,.4,.3,.2,.4};
    QList<double> X1 = X0;

    // quasimomentas
    QList<double> K = {-10.3,-2.24633, -0.714948, 0.714948, 2.24633,3.4,-1.2,2.6,6.1};

    // chain of position
    QList<QList<double> > Chain;


//*************************************************************************
//*************************************************************************



    ui->setupUi(this);
    //qDebug() << _sPermutator.getSpermSize(); // initially _N=0;
    _sPermutator.setSpermSize(N); // de facto particle number
    //qDebug() << _sPermutator.getSpermSize(); // N is established now!
    _sPermutator.sPermute(N,N); // all permutations preparation
    //_sPermutator.sProduct(N,0, c, X,K);
    //qDebug() << _sPermutator.SignFunction(1.2, 2.3);
    //_sPermutator.WaveFunction(N,c,X,K);
    //_sPermutator.Positions(X0,X1,delta,1,L);
    //_sPermutator.AddCollection(Chain, X1);
    //_sPermutator.Positions(X0,X1,delta,2,L);
    //_sPermutator.AddCollection(Chain, X1);
    _sPermutator.AddCollection(Chain, X0);
    qDebug() << Chain;
    for(int i = 1; i < N; i++)
    {
    _sPermutator.Metropolis(N,c,Chain,X0,X1,K,delta,N,L);
    qDebug() << X0 << i;
    }

    ui->tableWidget->setRowCount(Chain.size());
    ui->tableWidget->setColumnCount(N);

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
        double r = 255.*(Chain.at(i)).at(j);
        ui->tableWidget->setItem(i, j, new QTableWidgetItem);
        ui->tableWidget->item(i, j)->setBackground(QBrush(QColor(r ,255 - r ,r)));
        }
    }


}




MainWindow::~MainWindow()
{
    delete ui;

}
