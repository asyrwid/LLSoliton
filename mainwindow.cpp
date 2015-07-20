#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    qDebug() << _sPermutator.getSpermSize();
    _sPermutator.setSpermSize(5);
    _sPermutator.sPermute(4,4);


}

MainWindow::~MainWindow()
{
    delete ui;
}
