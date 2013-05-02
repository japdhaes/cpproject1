#ifndef JASTROW_H
#define JASTROW_H
#include <armadillo>
using namespace std;
using namespace arma;

class Jastrow
{
public:
    Jastrow();
    Jastrow(double _beta, int _nParticles);
    void setBeta(double _beta);
    void initialize(const mat &r);
    double evaluate(const mat &r);
    void setCurrentParticle(const int &i);
    double calcRatio();
    void acceptMove();
    void rejectMove();

    double getRatio();
    void setNewPosition(const mat &r);
protected:
    int nDimensions;
    int nParticles;
    int cp;

    double beta;
};

#endif // JASTROW_H
