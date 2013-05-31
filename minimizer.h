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
    void innercontraction();
    void outercontraction();

    void acceptReflection();
    void expansion();
    void decideCase();
    void calculateReflectedPoint();
    void calculateCenterOfGravity();

    void nelderMeadMethodDM();
    void reductionDM();
    void outercontractionDM();
    void innercontractionDM();
    void expansionDM();
    void calculateReflectedPointDM();
    void calculateCenterOfGravityDM();
    void decideCaseDM();
    void acceptReflectionDM();
protected:
    int myrank;
    int numprocs;
    int nParticles;

    double gammaNM;
    double rhoNM;
    double sigmaNM;
    double chiNM;

    double reflectedalpha, reflectedbeta, reflecteddist;
    double gravityalpha, gravitybeta, gravitydist;
    double expansionalpha, expansionbeta, expansiondist;
    double contractionalpha, contractionbeta, contractiondist;

    double energyReflectedPoint;
    double energyExpansionPoint;
    double energyContractionPoint;

    double totaldifference;
    double difference;

    mat points;
    VMCIS vmcsolver;
};

#endif // MINIMIZER_H
