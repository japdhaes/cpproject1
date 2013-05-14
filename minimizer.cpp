#include "minimizer.h"

Minimizer::Minimizer(int _myrank, int _numprocs, int nParticles):
    myrank(_myrank),
    numprocs(_numprocs),
    alphaNM(1.0),
    gammaNM(2.0),
    rhoNM(-0.5),
    sigmaNM(0.5),
    vmcsolver(VMCIS(_myrank, _numprocs, nParticles, 0, 0))
{
}

void Minimizer::bruteForce(int nParticles, double alphamin, double alphamax, double betamin, double betamax){
    double nrOfCycles=5e4;
    ofstream myfile;
    if(nParticles==2){
        myfile.open ("/home/jonathan/projectsCP/project1/data/helium/energies.txt");
    }
    else if(nParticles==4){
        myfile.open ("/home/jonathan/projectsCP/project1/data/beryllium/energies.txt");
    }
    else if(nParticles==10){
        myfile.open ("/home/jonathan/projectsCP/project1/data/neon/energies.txt");
    }
    else{
        cout << "Doing simulation with wrong number of particles. Exiting"<<endl;
        exit(1);
    }
    double alphastep=(alphamax-alphamin)/10;
    double betastep=(betamax-betamin)/10;
    for(double alpha=alphamin; alphamin<=alphamax; alpha+=alphastep){
        for(double beta=betamin; betamin<=betamax; beta+=betastep){
            VMCIS vmcsolver(myrank, numprocs, nParticles, alpha, beta);
            vmcsolver.setCycles(nrOfCycles);
            double energy=vmcsolver.runMonteCarloIntegration();
            myfile << alpha << " "<< beta << " "<<energy<<endl;
        }
    }
}

void Minimizer::nelderMeadMethod(){
    //points for the method, we start with a square + an initial good guess.
    //first column = alpha, second column = beta, third column is energy
    int nParticles=4;
    points.zeros(5,3);
//    VMCIS vmcsolver(myrank, numprocs, nParticles, 0, 0);
    vmcsolver.setCycles(1e5);
    double alpha, beta;
    //implemented for neon
    //point leftdown in square
    alpha=double(nParticles)/2; beta=0.0;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(0,0)=alpha;
    points(0,1)=beta;
    points(0,2)=vmcsolver.runMonteCarloIntegration();
    //point leftup in square
    alpha=double(nParticles)/2; beta=2.0;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(1,0)=alpha;
    points(1,1)=beta;
    points(1,2)=vmcsolver.runMonteCarloIntegration();
    //point rightup in square
    alpha=nParticles/2*3; beta=2.0;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(2,0)=alpha;
    points(2,1)=beta;
    points(2,2)=vmcsolver.runMonteCarloIntegration();
    //point rightdown in square
    alpha=nParticles/2*3; beta=0.0;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(3,0)=alpha;
    points(3,1)=beta;
    points(3,2)=vmcsolver.runMonteCarloIntegration();
    //point with good guess
    alpha=nParticles; beta=0.5;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(4,0)=alpha;
    points(4,1)=beta;
    points(4,2)=vmcsolver.runMonteCarloIntegration();

    //starting NM algorithm
    bool stop=false;
    while(!stop){

        //step 1 - order points
        orderPoints();
        cout << "POINTS "<<points;
        //step 2 - calculate center of gravity of all points except largest energy
        calculateCenterOfGravity();
        //step 3 - calculate reflection point
        calculateReflectedPoint();
        //step 4 - decide case dependent on the energy of the reflected point
        decideCase();
    }

}

//orders according to energy, ordering from lowest energy (row 0) to highest energy
void Minimizer::orderPoints(){
    for(int i=0; i<points.n_rows; i++){
        double minimumenergy=points(i,2);
        for(int j=i+1; j<points.n_rows;j++){
            if(points(j,2)<minimumenergy){
                minimumenergy=points(j,2);
            }
        }
        for(int j=i+1; j<points.n_rows;j++){
            if(points(j,2)==minimumenergy){
                rowvec temp = points.row(j);
                points.row(j)=points.row(i);
                points.row(i)=temp;
            }
        }
    }
}

void Minimizer::calculateCenterOfGravity(){
    gravityalpha=0.0;
    gravitybeta=0.0;
    for(int i=0; i<points.n_rows-1; i++){
        gravityalpha+=points(i,0);
        gravitybeta+=points(i,1);
    }
    gravityalpha/=(points.n_rows-1);
    gravitybeta/=(points.n_rows-1);
    cout << "gravityalpha "<<gravityalpha<< " gravitybeta "<<gravitybeta<<endl;
}

void Minimizer::calculateReflectedPoint(){
    reflectedalpha=(1.0+alphaNM)*gravityalpha-alphaNM*points(points.n_rows-1,0);
    reflectedbeta=(1.0+alphaNM)*gravitybeta-alphaNM*points(points.n_rows-1,1);
    if(reflectedbeta<0.0){
        reflectedbeta=0.0;
    }
    vmcsolver.setAlphaBeta(reflectedalpha,reflectedbeta);
    energyReflectedPoint=vmcsolver.runMonteCarloIntegration();
    cout << "reflectedalpha "<<reflectedalpha<<" beta "<<reflectedbeta<<" energy "<<energyReflectedPoint<<endl;
}

void Minimizer::decideCase(){
    bool found=false;
    int i=0;
    while(!found && i<points.n_rows){
        //if energy of i-th point is lower than the energy of the reflected point
        //we have not found the case so far
        if(points(i,2)<energyReflectedPoint){
            i++;
        }
        else{
            found=true;
        }
    }

    if(i==0){
        cout << "case expansion"<<endl;
        expansion();
    }
    else if(i<points.n_rows){
        cout << "case accept reflection"<<endl;
        acceptReflection();
    }
    else{
        cout << "case contraction"<<endl;
        contraction();
    }
}

void Minimizer::expansion(){
    expansionalpha=gravityalpha+gammaNM*(gravityalpha-points(points.n_rows-1,0));
    expansionbeta=gravitybeta+gammaNM*(gravitybeta-points(points.n_rows-1,1));
    if(expansionbeta<0.0){
        expansionbeta=0.0;
    }
    vmcsolver.setAlphaBeta(expansionalpha,expansionbeta);
    energyExpansionPoint=vmcsolver.runMonteCarloIntegration();
    cout << "expansionalpha= "<<expansionalpha<< " beta "<< expansionbeta << " energy "<<energyExpansionPoint<<endl;
    if(energyExpansionPoint<energyReflectedPoint){
        cout << "accepting expansionpoint"<<endl;
        points(points.n_rows-1,0)=expansionalpha;
        points(points.n_rows-1,1)=expansionbeta;
        points(points.n_rows-1,2)=energyExpansionPoint;
    }
    else{
        cout << "rejecting expansion point "<<endl;
        acceptReflection();
    }
}

void Minimizer::acceptReflection(){
    cout << "accepting reflection point "<<endl;
    cout << "reflectedalpha= "<<reflectedalpha<< " beta "<< reflectedbeta << " energy "<<energyReflectedPoint<<endl;
    points(points.n_rows-1,0)=reflectedalpha;
    points(points.n_rows-1,1)=reflectedbeta;
    vmcsolver.setAlphaBeta(points(points.n_rows-1,0),points(points.n_rows-1,1));
    points(points.n_rows-1,2)=energyReflectedPoint;
}

void Minimizer::contraction(){
    contractionalpha=gravityalpha+rhoNM*(gravityalpha-points(points.n_rows-1,0));
    contractionbeta=gravitybeta+rhoNM*(gravitybeta-points(points.n_rows-1,1));
    if(contractionbeta<0.0){
        contractionbeta=0.0;
    }
    vmcsolver.setAlphaBeta(contractionalpha,contractionbeta);
    energyContractionPoint=vmcsolver.runMonteCarloIntegration();

    //accepting contraction point
    if(energyContractionPoint<points(points.n_rows-1,2)){
        cout << "accepting contraction point "<<endl;
        cout << "contractionalpha= "<<contractionalpha<< " beta "<< contractionbeta << " energy "<<energyContractionPoint<<endl;
        points(points.n_rows-1,0)=contractionalpha;
        points(points.n_rows-1,1)=contractionbeta;
        points(points.n_rows-1,2)=energyContractionPoint;
    }
    else{
        reduction();
    }
}

void Minimizer::reduction(){
    cout << "reducting all points "<<endl;
    for(int i=1; i<points.n_rows;i++){
        points(i,0)=points(0,0)+sigmaNM*(points(i,0)-points(0,0));
        points(i,1)=points(0,1)+sigmaNM*(points(i,1)-points(0,1));
        if(points(i,1)<0){
            points(i,1)=0.0;
        }
        vmcsolver.setAlphaBeta(points(i,0),points(i,1));
        points(i,2)=vmcsolver.runMonteCarloIntegration();
    }
}
