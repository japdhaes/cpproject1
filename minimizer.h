#ifndef MINIMIZER_H
#define MINIMIZER_H
#include "vmcis.h"
#include "vmcbf.h"
class Minimizer
{
public:
    void bruteForce(double alphamin, double alphamax, double betamin, double betamax);

    Minimizer(int _myrank, int _numprocs, int nParticles);
    void bruteForce(int nParticles, double alphamin, double alphamax, double betamin, double betamax);
    void nelderMeadMethod();
    void orderPoints();
    void reduction();
    void contraction();

    void acceptReflection();
    void expansion();
    void decideCase();
    void calculateReflectedPoint();
    void calculateCenterOfGravity();
protected:
    int myrank;
    int numprocs;

    double alphaNM;
    double gammaNM;
    double rhoNM;
    double sigmaNM;

    double reflectedalpha, reflectedbeta;
    double gravityalpha, gravitybeta;
    double expansionalpha, expansionbeta;
    double contractionalpha, contractionbeta;

    double energyReflectedPoint;
    double energyExpansionPoint;
    double energyContractionPoint;

    double totaldifference;
    double difference;

    mat points;
    VMCIS vmcsolver;
};

#endif // MINIMIZER_H
