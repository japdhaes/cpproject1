TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    vmcsolver.cpp \
    lib.cpp \
    mainapplication.cpp \
    slater.cpp \
    orbital.cpp \
    jastrow.cpp \
    vmcis.cpp \
    vmcbf.cpp \
    testdirectory.cpp \
    datalogger.cpp \
    minimizer.cpp \
    dimoleculewavefunction.cpp \
    dimoleculeorbitals.cpp \
    orbitals.cpp \
    atomwavefunction.cpp

HEADERS += \
    vmcsolver.h \
    lib.h \
    mainapplication.h \
    slater.h \
    orbital.h \
    jastrow.h \
    vmcis.h \
    vmcbf.h \
    testdirectory.h \
    datalogger.h \
    minimizer.h \
    dimoleculewavefunction.h \
    dimoleculeorbitals.h \
    orbitals.h \
    atomwavefunction.h

LIBS+= -larmadillo

QMAKE_CXXFLAGS_WARN_ON += -Wno-reorder

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}
