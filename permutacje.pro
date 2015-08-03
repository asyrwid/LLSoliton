#-------------------------------------------------
#
# Project created by QtCreator 2015-07-20T21:50:51
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

TARGET = permutacje
TEMPLATE = app
CONFIG += c++11

SOURCES += main.cpp\
        mainwindow.cpp \
    spermutator.cpp \
    permutator.cpp

HEADERS  += mainwindow.h \
    spermutator.h \
    permutator.h

FORMS    += mainwindow.ui

QMAKE_CXXFLAGS  += -O3
