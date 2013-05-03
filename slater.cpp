#include "slater.h"

//this my trial wave function class
//Slater::Slater(int myrank, int numprocs, double _alpha, double _beta):
//    VMCSolver(myrank, numprocs)
//{
//    this->charge=4;
//    this->nParticles=this->charge;
//    this->alpha=_alpha;
//    this->beta=_beta;
//    orbs=Orbital(this->alpha);
//}

Slater::Slater(double _alpha, int _nParticles){
    this->alpha=_alpha;
    this->nDimensions=3;
    this->nParticles=_nParticles;
    this->hp=_nParticles/2;
    this->orbs=Orbital(this->alpha);
    this->rOld = zeros<mat>(nParticles, nDimensions);
    this->rNew = zeros<mat>(nParticles, nDimensions);
    this->sdup = zeros<mat>(hp, hp);
    this->sddown = zeros<mat>(hp, hp);
    this->sdupnew = zeros<mat>(hp, hp);
    this->sddownnew = zeros<mat>(hp, hp);
    this->sdupinverse = zeros<mat>(hp, hp);
    this->sddowninverse = zeros<mat>(hp, hp);
    this->sdupinversenew = zeros<mat>(hp, hp);
    this->sddowninversenew = zeros<mat>(hp, hp);
    this->R=1.0;
    this->S=zeros(1,this->hp);
}

//debugged
void Slater::initialize(const mat &r){
    rNew = rOld = r;
    //===============================
    //initializing inverse matrices
    //===============================
    vec3 temp;
    for(int particle=0; particle<this->hp; particle++){
        for(int orbital=0; orbital<this->hp;orbital++){
            temp = rNew.row(particle);
            this->sdup(particle,orbital)=orbs.hydrogenWF(temp, orbital);
            temp = rNew.row(particle+this->hp);
            this->sddown(particle,orbital)=orbs.hydrogenWF(temp, orbital);
        }
    }
    this->sdupinverse=inv(sdup);
    this->sddowninverse=inv(sddown);

    this->sdupnew = this->sdup;
    this->sddownnew = this->sddown;
    this->sdupinversenew = this->sdupinverse;
    this->sddowninversenew = this->sddowninverse;
}

double Slater::localEnergyClosedForm(const mat &r)
{
    cout << __LINE__<<endl;
    //impossible
    return 0;
}

double Slater::evaluate(){
    return det(sdup)*det(sddown);
}

//slaterdeterminant
//debugged
double Slater::evaluate(const mat &r){
    mat sduptemp = zeros<mat>(hp, hp);
    mat sddowntemp = zeros<mat>(hp, hp);
    vec3 temp;
    for(int i=0; i<hp; i++){
        for(int j=0; j<hp;j++){
            temp = r.row(j);
            sduptemp(i,j)=orbs.hydrogenWF(temp, i);
            temp = r.row(j+hp);
            sddowntemp(i,j)=orbs.hydrogenWF(temp, i);
        }
    }
    return det(sduptemp)*det(sddowntemp);
}

double Slater::getRatio(){
    return this->R;
}

//debugged
double Slater::calcRatio(){
    //formula 68 on slides
    //R = sum_j d_ij(rNEW!) d_ji ^-1 (rOLD!)
    vec3 currentpos = rNew.row(this->cp);
    this->R=0.0;

    //spin UP
    if(spinup){
        int i=this->cp;
        for(int j=0; j < this->hp; j++){
            this->R+=this->sdupnew(i,j)*this->sdupinverse(j,i);
        }
//        cout << "R= "<<this->R << endl << "sdupnew "<<this->sdupnew<<"sdupinverseold"<< this->sdupinverse<<endl;
    }
    //spin DOWN
    else{
        int i=this->cp-this->hp;
        for(int j=0; j < this->hp; j++){
            this->R+=this->sddownnew(i,j)*this->sddowninverse(j,i);
        }
    }

    return this->R;
}

//debugged
void Slater::acceptMove(){
    this->rOld.row(this->cp)=this->rNew.row(this->cp);
    this->setNewSlaterInverse();
    this->updateSlaterAndInverse();
}

void Slater::setAlpha(const double _alpha)
{
    this->alpha=_alpha;
}

//debugged
void Slater::setCurrentParticle(const int &i)
{
    this->cp=i;
    //spinUP
    if(this->cp<this->hp){
        this->spinup=true;
    }
    //spinDOWN
    else{
        this->spinup=false;
    }
}

//debugged
void Slater::setNewPosition(const mat &r){
    this->rNew=r;
    this->setSlaterNew();
}

//debugged
void Slater::setSlaterNew(){
    //spin UP
    if(spinup){
        int i = this->cp;
        for(int j=0; j<this->hp;j++){
            rowvec temp = rNew.row(i);
            this->sdupnew(i,j) = orbs.hydrogenWF(temp,j);
        }
    }
    //spin DOWN
    else {
        int i = this->cp-this->hp;
        for(int j=0; j<this->hp;j++){
            rowvec temp = rNew.row(this->cp);
            this->sddownnew(i,j) = orbs.hydrogenWF(temp,j);
        }
    }
}

//debugged
void Slater::updateSlaterAndInverse(){
    //spin UP
    if(spinup){
        int i=this->cp;
        this->sdup.row(i) = this->sdupnew.row(i);
        this->sdupinverse = this->sdupinversenew;
    }
    //spin DOWN
    else{
        int i=this->cp-this->hp;
        this->sddown.row(i) = this->sddownnew.row(i);
        this->sddowninverse = this->sddowninversenew;
    }
}

//debugged
void Slater::rejectMove(){
    this->rNew.row(this->cp) = this->rOld.row(this->cp);
    //spin UP
    if(spinup){
        int i = this->cp;
        sdupnew.row(i) = sdup.row(i);
    }
    //spin DOWN
    else{
        int i=this->cp-this->hp;
        this->sddownnew.row(i) = this->sddown.row(i);
    }
}

//debugged
void Slater::setNewSlaterInverse(){
    //ratio should already be calculated when checking to accept or reject the move!
    this->R=this->getRatio();
    this->S=zeros(1,this->hp);

    //spin UP
    if(this->cp<this->hp){
        int i= this->cp;
        //loop over all other up-particles
        for(int j=0; j<this->hp; j++){
            for(int k=0; k<this->hp; k++){
                //formula 69 on slides
                //     d_ik(rNEW!)         d_kj^-1(rOLD!)
                S(j)+=sdupnew(i,k)*this->sdupinverse(k,j);
            }
        }
        for(int j=0; j<this->hp;j++){
            for(int k=0; k<this->hp;k++){
                if(k!=i){
                    //formula 70 on slides
                    //                  d_kj^-1(rOLD!)                d_ki^-1(rOLD!)
                    sdupinversenew(j,k)=sdupinverse(j,k) - (S(k)/R)*sdupinverse(j,i);
                }
                else{
                    //formula 71 on slides
                    //                          d_ki^-1(rOLD!)
                    sdupinversenew(j,k)=(1.0/R)*sdupinverse(j,i);
                }
            }
        }
    }
    //spin DOWN
    else{
        int i=this->cp-this->hp;
        //loop over all other down-particles
        for(int j=0; j<this->hp; j++){
            for(int k=0; k<this->hp; k++){
                //formula 69 on slides
                //     d_ik(rNEW!)         d_kj^-1(rOLD!)
                S(j)+=sddownnew(i,k)*this->sddowninverse(k,j);
            }
        }
        for(int j=0; j<this->hp;j++){
            for(int k=0; k<this->hp;k++){
                if(k!=i){
                    //formula 70 on slides
                    //                  d_kj^-1(rOLD!)                d_ki^-1(rOLD!)
                    sddowninversenew(j,k)=sddowninverse(j,k) - (S(k)/R)*sddowninverse(j,i);
                }
                else{
                    //formula 71 on slides
                    //                          d_ki^-1(rOLD!)
                    sddowninversenew(j,k)=(1.0/R)*sddowninverse(j,i);
                }
            }
        }
    }

}

