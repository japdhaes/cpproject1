#include "slater.h"

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
    this->rOld=this->rNew;
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

void Slater::setAlpha(const double _alpha)
{
    this->alpha=_alpha;
    this->orbs.setAlpha(_alpha);
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
    this->setNewSlaterInverse();
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
void Slater::rejectMove(){
    this->rNew.row(this->cp) = this->rOld.row(this->cp);
    //spin UP
    if(spinup){
        int i = cp;
        this->sdupnew.row(i) = sdup.row(i);
        this->sdupinversenew=this->sdupinverse;
    }
    //spin DOWN
    else{
        int i=this->cp-this->hp;
        this->sddownnew.row(i) = this->sddown.row(i);
        this->sddowninversenew=this->sddowninverse;
    }
}

//debugged
void Slater::setNewSlaterInverse(){
    this->R=this->calcRatio();
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

rowvec Slater::localGradient(const int &i){
    //formula 72 on slides
    rowvec answer(this->nDimensions);
    answer.zeros();
    rowvec temp=this->rNew.row(i);
    bool spinup=(i<this->hp);
    if(spinup){
        for(int j=0; j<this->hp; j++){
            //      gradient_i phi_j(r_i) *           d_ji^-1 (rNEW)
            answer+=orbs.gradient(temp,j)*this->sdupinversenew(j,i);
        }
    }
    else{
        for(int j=0; j<this->hp; j++){
            answer+=orbs.gradient(temp,j)*this->sddowninversenew(j,i-this->hp);
        }
    }
    return answer;
}

double Slater::localLaplacian(const int &i){
    //formula 73 on slides
    bool spinup=(i<this->hp);
    double answer=0.0;
    rowvec temp=this->rNew.row(i);
    if(spinup){
        for(int j=0; j<this->hp; j++){
            answer+=orbs.laplacian(temp,j)*this->sdupinversenew(j,i);
        }
    }
    else{
        for(int j=0; j<this->hp; j++){
            answer+=orbs.laplacian(temp,j)*this->sddowninversenew(j,i-this->hp);
        }
    }
    return answer;
}

mat Slater::gradient(const mat &r, const double &h)
{
    mat rPlus, rMinus;
    rPlus = rMinus = r;
    mat temp(nParticles, nDimensions);
    double wfCurrent = evaluate(r);
    double dfactor = 1.0/(wfCurrent*2.0*h);
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            double wfMinus = evaluate(rMinus);
            double wfPlus = evaluate(rPlus);
            temp(i,j) = (wfPlus - wfMinus)*dfactor;
            rPlus(i,j) = rMinus(i,j) = r(i,j);
        }
    }

    return temp;
}

double Slater::alphaGradient(const int &i)
{
    /* Calculates the gradient for particle i */
    double delta = 0.0;
    if (i < this->hp){
        for (int j = 0; j < this->hp; j++){
            delta += orbs.alphaGradient(rNew.row(i),j)*this->sdupinversenew(j,i);
        }
    }
    else{
        for (int j = 0; j < this->hp; j++){
            delta += orbs.alphaGradient(rNew.row(i),j)*this->sddowninversenew(j,i-this->hp);
        }
    }

    return delta;
}
