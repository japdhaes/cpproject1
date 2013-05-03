#include "jastrow.h"

Jastrow::Jastrow(double _beta, int _nParticles):
    a(zeros<mat>(_nParticles,_nParticles)),
    beta(_beta),
    nParticles(_nParticles),
    nDimensions(3),
    rOld(zeros<mat>(nParticles, nDimensions)),
    rNew(zeros<mat>(nParticles, nDimensions)),
    rij(zeros<mat>(_nParticles,_nParticles)),
    rijNew(zeros<mat>(_nParticles,_nParticles)),
    fij(zeros<mat>(_nParticles,_nParticles)),
    fijNew(zeros<mat>(_nParticles,_nParticles))
{
    initialize_a();
}

void Jastrow::setBeta(double _beta)
{
    this->beta=_beta;
}

void Jastrow::initialize_a(){
    for(int i=0; i<this->nParticles;i++){
        for(int j=0; j<this->nParticles;j++){
            if((i+j)%2){
                a(i,j)=0.5;
            }
            else{
                a(i,j)=0.25;
            }
        }
    }
}

//debugged
void Jastrow::initialize(const mat &r)
{
    rNew = rOld = r;
    rijNew=rij=calc_rij(rNew);
    fijNew=fij=calc_fij(rijNew);
}

mat Jastrow::calc_rij(const mat &r){
    mat answer(r.n_rows,r.n_rows);
    answer.zeros();
    for(int i=0; i<this->nParticles;i++){
        for(int j=0; j<i; j++){
            for(int k=0; k<this->nDimensions;k++){
                answer(i,j) += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            answer(i,j)=sqrt(answer(i,j));
            answer(j,i)=answer(i,j);
        }
    }
    return answer;
}

mat Jastrow::calc_fij(const mat &_rij){
    mat answer(_rij);
    double temp;
    for(int i=0; i<this->nParticles;i++){
        for(int j=0;j<this->nParticles;j++){
            temp=_rij(i,j);
            answer(i,j) = a(i,j)*temp/(1.0+this->beta*temp);
            answer(j,i) = a(i,j)*temp/(1.0+this->beta*temp);
        }
    }
    return answer;
}

//void Jastrow::calculate_fij()
//{
//    double rij;
//    for (int i = 0; i < nParticles; i++)
//    {
//        for (int j = 0; j < nParticles; j++)
//        {
//            rij = rijOld(i,j);
//            fijOld(i,j) = a(i,j)*rij/(1.0 + beta*rij);
//        }
//    }
//}

double Jastrow::evaluate(const mat &r)
{
    mat rijmatrix=this->calc_rij(r);
    mat fijmatrix=this->calc_fij(rijmatrix);
    double arg=0.0;
    for(int i=0; i<this->nParticles;i++){
        for(int j=0; j<i; j++){
            arg+=fijmatrix(i,j);
        }
    }
    return exp(arg);
}

void Jastrow::setCurrentParticle(const int &i)
{
    this->cp=i;
}

void Jastrow::setNewPosition(const mat &r){
    this->rNew=r;
    this->setNewRij();
    this->setNewFij();

}

void Jastrow::setNewRij(){
    int i=this->cp;
//    cout << rNew<<endl;
    for(int j=0; j<i; j++){
        rijNew(i,j)=0.0;
        for(int k=0; k<this->nDimensions;k++){
            rijNew(i,j)+=(rNew(i,k)-rNew(j,k))*(rNew(i,k)-rNew(j,k));
        }
        rijNew(i,j)=sqrt(rijNew(i,j));
        rijNew(j,i)=rijNew(i,j);
    }
    for(int j=i+1; j<this->nParticles; j++){
        rijNew(i,j)=0.0;
        for(int k=0; k<this->nDimensions;k++){
            rijNew(i,j)+=(rNew(i,k)-rNew(j,k))*(rNew(i,k)-rNew(j,k));
        }
        rijNew(i,j)=sqrt(rijNew(i,j));
        rijNew(j,i)=rijNew(i,j);
    }
}

void Jastrow::setNewFij(){
    int i=this->cp;
    double temp;
    for(int j=0; j<i;j++){
        temp=rijNew(i,j);
        fijNew(i,j)=a(i,j)*temp/(1.0+this->beta*temp);
        fijNew(j,i)=fijNew(i,j);
    }
    for(int j=i+1; j<this->nParticles;j++){
        temp=rijNew(i,j);
        fijNew(i,j)=a(i,j)*temp/(1.0+this->beta*temp);
        fijNew(j,i)=fijNew(i,j);
    }
}

double Jastrow::calcRatio()
{
    double arg=0.0;
    int i=this->cp;
//    cout << "new"<<fijNew<<endl;
//    cout <<"old"<<fij<<endl;
    for(int j=0; j<i; j++){
        arg+=fijNew(i,j)-fij(i,j);
    }
    for(int j=i+1; j<this->nParticles;j++){
        arg+=fijNew(i,j)-fij(i,j);
    }
    this->R=exp(arg);
//    cout << this->R*this->R <<endl;
    return this->R;
}


double Jastrow::getRatio(){
    return this->R;
}

void Jastrow::acceptMove()
{
    rOld.row(this->cp)=rNew.row(this->cp);
    rij.row(this->cp)=rijNew.row(this->cp);
    rij.col(this->cp)=rijNew.col(this->cp);
    fij.row(this->cp)=fijNew.row(this->cp);
    fij.col(this->cp)=fijNew.col(this->cp);
}


void Jastrow::rejectMove()
{
    rNew.row(this->cp)=rOld.row(this->cp);
    rijNew.row(this->cp)=rij.row(this->cp);
    rijNew.col(this->cp)=rij.col(this->cp);
    fijNew.row(this->cp)=fij.row(this->cp);
    fijNew.col(this->cp)=fij.col(this->cp);
}
