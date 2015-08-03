#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QCloseEvent>
#include <QTableWidget>
#include <QElapsedTimer>
#include <iostream>
#include <fstream>
/* srand example */
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <omp.h>       /* OpenMP */




MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

    QElapsedTimer timer;
    timer.start();


    srand (time(NULL));
    qDebug() << "First number: %d\n" << rand()%100;
// initial values!
//*************************************************************************
//*************************************************************************


    int    N = 8; // number of particles
    double c = 24.; // coupling constant
    double L = 1; // system length
    double delta = 0.3; // maximal "jump" between single particle
                        // position from one realization to another

    // initial collection of particle positions
    QList<double> X0 = {.1,.12,.13,.122,.132,.111,.113,.118};
    QList<double> X1 = X0;

    // quasimomentas
    QList<double> K = {-0.752974,-0.30724,0.107856,0.552964,5.73022,6.17533,6.59043,7.03616};

    // chain of position
    QList<QList<double> > Chain;

    //chain of wave function values
    QList<double> WaveFValues;


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
    for(int i = 0; i < 10000; i++)
    {
    _sPermutator.Metropolis(N,c,Chain,WaveFValues,X0,X1,K,delta,N,L);
    qDebug() << i;
    //qDebug() << WaveFValues;
    }
    for(int i = 0; i < Chain.size(); i++)
    {
    qDebug() << Chain.at(i);
    }

    qDebug() << Chain ;




    std::ofstream create ("test.txt");
    std::fstream fs;
    fs.open ("test.txt");
    for(int i = 0; i < Chain.size(); i++)
    {
        for(int j = 0; j < N; j++ )
        {
            fs << (Chain.at(i)).at(j) << " ";
        }
        fs << "\n";
    }
    fs.close();


    ui->tableWidget->setRowCount(Chain.size());
    ui->tableWidget->setColumnCount(N);

    for(int i = 0; i < Chain.size(); i++)
    {
        for(int j = 0; j < N; j++)
        {
        double r = 16.*(Chain.at(i)).at(j);
        ui->tableWidget->setItem(i, j, new QTableWidgetItem);
        ui->tableWidget->item(i, j)->setBackground(QBrush(QColor(0  ,255 - r*r ,r*r)));
        }
    }
qDebug() << (timer.elapsed())/1000.;

}




MainWindow::~MainWindow()
{
    delete ui;

}
