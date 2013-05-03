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
    void initialize_a();
    mat calc_rij(const mat &r);
    mat calc_fij(const mat &_rij);
    void setNewRij();
    void setNewFij();
protected:
    int nDimensions;
    int nParticles;
    int cp;

    double beta;
    double R;

    mat a;

    mat rij;
    mat rijNew;
    mat fij;
    mat fijNew;
    mat rOld;
    mat rNew;
};

#endif // JASTROW_H