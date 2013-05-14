#ifndef SLATER_H
#define SLATER_H
#include "orbital.h"
#include <armadillo>
using namespace std;
using namespace arma;

class Slater
{
public:
    Slater();
//    Slater(int myrank, int numprocs, double _alpha, double _beta);
    Slater(double _alpha, int _nParticles);
    void setAlpha(const double _alpha);
    void setCurrentParticle(const int &i);
    void setNewPosition(const mat &r);
    void setSlaterNew();
    void setNewSlaterInverse();

    double localEnergyClosedForm(const mat &r);

    double evaluate();
    double evaluate(const mat &r);

    double getRatio();
    double calcRatio();

    void initialize(const mat &r);
    void acceptMove();
    void rejectMove();
    void updateSlaterAndInverse();


    rowvec localGradient(const int &i);
    double localLaplacian(const int &i);
    double localLaplacian2(const int &i);
    mat gradient(const mat &r, const double &h);
    double alphaGradient(const int &i);
protected:
    int nDimensions;
    int nParticles;
    //half particles
    int hp;
    //current particle
    int cp;
    //spinup? true=spinup, false=spindown
    bool spinup;

    mat rNew, rOld;
    mat sdup, sddown;
    mat sdupinverse, sddowninverse;
    mat sdupnew, sddownnew;
    mat sdupinversenew, sddowninversenew;
    Orbital orbs;

    double alpha;
    double R;
    vec3 S;
};

#endif // SLATER_H
