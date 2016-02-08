QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    permutator.cpp \
    decomposition.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    permutator.h \
    decomposition.h \
    utils.h

QMAKE_CXXFLAGS  += -O3
