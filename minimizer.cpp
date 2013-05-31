#include "minimizer.h"

Minimizer::Minimizer(int _myrank, int _numprocs, int nParticles):
    myrank(_myrank),
    numprocs(_numprocs),
    rhoNM(1.0),
    chiNM(2.0),
    gammaNM(0.5),
    sigmaNM(0.5),
    vmcsolver(VMCIS(_myrank, _numprocs, nParticles, 0, 0)),
    nParticles(nParticles)
{
}

void Minimizer::bruteForce(double alphamin, double alphamax, double betamin, double betamax){
    double nrOfCycles=1e5;
    cout << "starting bruteforce minimizer with "<<nParticles <<" particles, alphamin=" <<alphamin;
    cout <<" , alphamax="<<alphamax << " , betamin="<<betamin << " , betamax="<<betamax<<endl;
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
    double alphastep=(alphamax-alphamin)/20;
    double betastep=(betamax-betamin)/15;
    for(double alpha=alphamin; alpha<=alphamax; alpha+=alphastep){
        for(double beta=betamin; beta<=betamax; beta+=betastep){
            VMCIS vmcsolver(myrank, numprocs, nParticles, alpha, beta);
            cout << "simulating alpha="<<alpha<< " beta="<<beta<<endl;
            vmcsolver.setCycles(nrOfCycles);
            double energy=vmcsolver.runMonteCarloIntegration();
            myfile << alpha << " "<< beta << " "<<energy<<endl;
        }
        myfile << endl;
    }
}

void Minimizer::nelderMeadMethod(){
    double tolerance=1e-6;
    //points for the method, we start with a square + an initial good guess.
    //first column = alpha, second column = beta, third column is energy
    //int nParticles=this->;
    points.zeros(3,3);
//    VMCIS vmcsolver(myrank, numprocs, nParticles, 0, 0);
    vmcsolver.setCycles(1e6);
    double alpha, beta;
    //implemented for neon
    //point leftdown in square
    alpha=nParticles; beta=0.5;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(0,0)=alpha;
    points(0,1)=beta;
    points(0,2)=vmcsolver.runMonteCarloIntegration();
//    //point leftup in square
//    alpha=double(nParticles)/2; beta=2.0;
//    vmcsolver.setAlphaBeta(alpha, beta);
//    points(1,0)=alpha;
//    points(1,1)=beta;
//    points(1,2)=vmcsolver.runMonteCarloIntegration();
    //point rightup in square
    alpha=nParticles*0.75; beta=0.25;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(1,0)=alpha;
    points(1,1)=beta;
    points(1,2)=vmcsolver.runMonteCarloIntegration();
//    //point rightdown in square
//    alpha=nParticles/2*3; beta=0.0;
//    vmcsolver.setAlphaBeta(alpha, beta);
//    points(3,0)=alpha;
//    points(3,1)=beta;
//    points(3,2)=vmcsolver.runMonteCarloIntegration();
    //point with good guess
    alpha=nParticles*0.9; beta=0.15;
    vmcsolver.setAlphaBeta(alpha, beta);
    points(2,0)=alpha;
    points(2,1)=beta;
    points(2,2)=vmcsolver.runMonteCarloIntegration();

    //starting NM algorithm
    bool stop=false;
    double averageEnergy=0.0;
    double previousEnergy=0.0;
    double previousAverageAlpha=0.0;
    double previousAverageBeta=0.0;
    double currentAverageAlpha=0.0;
    double currentAverageBeta=0.0;
    difference=1.0;
    while(difference>tolerance){
        //step 1 - order points
        orderPoints();
        if(myrank==0){
            cout.precision(11);
            cout.setf(ios::fixed);

            points.raw_print(cout, "POINTS=");
        //cout << "POINTS "<<points;
        }

        //step 2 - calculate center of gravity of all points except largest energy
        calculateCenterOfGravity();
        //step 3 - calculate reflection point
        calculateReflectedPoint();
        //step 4 - decide case dependent on the energy of the reflected point
        decideCase();

        //deciding to stop the algorithm
        averageEnergy=0.0;
        for(int i=0; i<points.n_rows;i++){
            currentAverageAlpha+=points(i,0);
            currentAverageBeta+=points(i,1);
        }
        currentAverageAlpha/=points.n_rows;
        currentAverageBeta/=points.n_rows;

        difference=sqrt(pow((currentAverageAlpha-previousAverageAlpha),2)+pow((currentAverageBeta-previousAverageBeta),2));

//        if(difference<tolerance){
//            cout << "I am rank "<<myrank<<" and I reached tolerance level "<< difference<<endl;
//        }
        previousAverageAlpha=currentAverageAlpha;
        previousAverageBeta=currentAverageBeta;
        cout << "finishing iteration "<<endl;
    }
    cout << "out of loop "<<endl;
}

//orders according to energy, ordering from lowest energy (row 0) to highest energy
void Minimizer::orderPoints(){
    for(int i=0; i<points.n_rows; i++){
        double minimumenergy=points(i,points.n_cols-1);
        for(int j=i+1; j<points.n_rows;j++){
            if(points(j,points.n_cols-1)<minimumenergy){
                minimumenergy=points(j,points.n_cols-1);
            }
        }
        for(int j=i+1; j<points.n_rows;j++){
            if(points(j,points.n_cols-1)==minimumenergy){
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
}

void Minimizer::calculateReflectedPoint(){
    reflectedalpha=(1.0+rhoNM)*gravityalpha-rhoNM*points(points.n_rows-1,0);
    reflectedbeta=(1.0+rhoNM)*gravitybeta-rhoNM*points(points.n_rows-1,1);

    if(reflectedbeta<0.0){
        reflectedbeta=0.39;
    }
    vmcsolver.setAlphaBeta(reflectedalpha,reflectedbeta);
    energyReflectedPoint=vmcsolver.runMonteCarloIntegration();
}

void Minimizer::decideCase(){
    bool found=false;
    int i=0;
    while(!found ){
        //if energy of i-th point is lower than the energy of the reflected point
        //we have not found the case so far
        if(points(i,points.n_cols-1)<energyReflectedPoint){
            i++;
        }
        else if (i==points.n_rows){
            i++;
            found=true;
        }
        else{
            found=true;
        }
    }

    if(i==0){
        expansion();
    }
    else if(i<points.n_rows-1){
        acceptReflection();
    }
    else if(i==points.n_rows-1){
        outercontraction();
    }
    else{
        innercontraction();
    }
}
void Minimizer::decideCaseDM(){
    bool found=false;
    int i=0;
    while(!found ){
        //if energy of i-th point is lower than the energy of the reflected point
        //we have not found the case so far
        if(points(i,points.n_cols-1)<energyReflectedPoint){
            i++;
        }
        else if (i==points.n_rows){
            i++;
            found=true;
        }
        else{
            found=true;
        }
    }

    if(i==0){
        expansionDM();
    }
    else if(i<points.n_rows-1){
        acceptReflectionDM();
    }
    else if(i==points.n_rows-1){
        outercontractionDM();
    }
    else{
        innercontractionDM();
    }
}

void Minimizer::expansion(){
    expansionalpha=gravityalpha-chiNM*(gravityalpha-reflectedalpha);
    expansionbeta=gravitybeta-chiNM*(gravitybeta-reflectedbeta);
    if(expansionbeta<0.0){
        expansionbeta=0.39;
    }
    vmcsolver.setAlphaBeta(expansionalpha,expansionbeta);
    energyExpansionPoint=vmcsolver.runMonteCarloIntegration();

    if(energyExpansionPoint<energyReflectedPoint){
        if(myrank==0)
        cout <<"accepting expansion "<<endl;
        points(points.n_rows-1,0)=expansionalpha;
        points(points.n_rows-1,1)=expansionbeta;
        points(points.n_rows-1,2)=energyExpansionPoint;
    }
    else{
        acceptReflection();
    }
}


void Minimizer::acceptReflection(){
    if(myrank==0)
    cout << "accept reflection "<<endl;
    points(points.n_rows-1,0)=reflectedalpha;
    points(points.n_rows-1,1)=reflectedbeta;
    vmcsolver.setAlphaBeta(points(points.n_rows-1,0),points(points.n_rows-1,1));
    points(points.n_rows-1,2)=energyReflectedPoint;
}

void Minimizer::innercontraction(){

    contractionalpha=gravityalpha-gammaNM*(gravityalpha-points(points.n_rows-1,0));
    contractionbeta=gravitybeta-gammaNM*(gravitybeta-points(points.n_rows-1,1));
    if(contractionbeta<0.0){
        contractionbeta=0.39;
    }
    vmcsolver.setAlphaBeta(contractionalpha,contractionbeta);
    energyContractionPoint=vmcsolver.runMonteCarloIntegration();

    //accepting contraction point
    if(energyContractionPoint<points(points.n_rows-1,2)){
        if(myrank==0)
        cout << "accept innercontraction "<<endl;
        points(points.n_rows-1,0)=contractionalpha;
        points(points.n_rows-1,1)=contractionbeta;
        points(points.n_rows-1,2)=energyContractionPoint;
    }
    else{
        reduction();
    }
}

void Minimizer::outercontraction()
{
    contractionalpha=gravityalpha+gammaNM*(reflectedalpha-gravityalpha);
    contractionbeta=gravitybeta+gammaNM*(reflectedbeta-gravitybeta);
    if(contractionbeta<0.0){
        contractionbeta=0.39;
    }
    vmcsolver.setAlphaBeta(contractionalpha,contractionbeta);
    energyContractionPoint=vmcsolver.runMonteCarloIntegration();

    //accepting contraction point
    if(energyContractionPoint<points(points.n_rows-1,2)){
        if(myrank==0)
        cout << "accept outercontraction "<<endl;
        points(points.n_rows-1,0)=contractionalpha;
        points(points.n_rows-1,1)=contractionbeta;
        points(points.n_rows-1,2)=energyContractionPoint;
    }
    else{
        reduction();
    }
}

void Minimizer::reduction(){
    if(myrank==0)
    cout << "reduction "<<endl;
    for(int i=1; i<points.n_rows;i++){
        points(i,0)=points(0,0)+sigmaNM*(points(i,0)-points(0,0));
        points(i,1)=points(0,1)+sigmaNM*(points(i,1)-points(0,1));
        if(points(i,1)<0){
            points(i,1)=0.0;
        }
        vmcsolver.setAlphaBeta(points(i,0),points(i,1));
        double reductedEnergy=vmcsolver.runMonteCarloIntegration();
        points(i,2)=reductedEnergy;
    }
}

void Minimizer::nelderMeadMethodDM(){
    double tolerance=1e-6;
    //points for the method, we start with a square + an initial good guess.
    //first column = alpha, second column = beta, third column is energy
    //int nParticles=this->;
    points.zeros(4,4);
//    VMCIS vmcsolver(myrank, numprocs, nParticles, 0, 0);
    vmcsolver.setCycles(1e6);
    double alpha, beta, dist;
    //implemented for neon
    //point leftdown in square
    alpha=nParticles; beta=0.5; dist=0.5;
    vmcsolver.setAlphaBeta(alpha, beta);
    vmcsolver.setdist(dist);
    points(0,0)=alpha;
    points(0,1)=beta;
    points(0,2)=dist;
    points(0,3)=vmcsolver.runMonteCarloIntegration();
//    //point leftup in square
//    alpha=double(nParticles)/2; beta=2.0;
//    vmcsolver.setAlphaBeta(alpha, beta);
//    points(1,0)=alpha;
//    points(1,1)=beta;
//    points(1,2)=vmcsolver.runMonteCarloIntegration();
    //point rightup in square
    alpha=nParticles*0.75; beta=0.25; dist=1.0;
    vmcsolver.setAlphaBeta(alpha, beta);
    vmcsolver.setdist(dist);
    points(1,0)=alpha;
    points(1,1)=beta;
    points(1,2)=dist;
    points(1,3)=vmcsolver.runMonteCarloIntegration();
//    //point rightdown in square
//    alpha=nParticles/2*3; beta=0.0;
//    vmcsolver.setAlphaBeta(alpha, beta);
//    points(3,0)=alpha;
//    points(3,1)=beta;
//    points(3,2)=vmcsolver.runMonteCarloIntegration();
    //point with good guess
    alpha=nParticles*0.9; beta=0.15; dist=1.5;
    vmcsolver.setAlphaBeta(alpha, beta);
    vmcsolver.setdist(dist);
    points(2,0)=alpha;
    points(2,1)=beta;
    points(2,2)=dist;
    points(2,3)=vmcsolver.runMonteCarloIntegration();

    alpha=nParticles*1.1; beta=0.4; dist=1.5;
    vmcsolver.setAlphaBeta(alpha, beta);
    vmcsolver.setdist(dist);
    points(3,0)=alpha;
    points(3,1)=beta;
    points(3,2)=dist;
    points(3,3)=vmcsolver.runMonteCarloIntegration();



    //starting NM algorithm
    bool stop=false;
    double averageEnergy=0.0;
    double previousEnergy=0.0;
    double previousAverageAlpha=0.0;
    double previousAverageBeta=0.0;
    double previousAverageDist=0.0;
    double currentAverageAlpha=0.0;
    double currentAverageBeta=0.0;
    double currentAverageDist=0.0;
    difference=1.0;
    while(difference>tolerance){
        //step 1 - order points
        orderPoints();
        if(myrank==0){
            cout.precision(11);
            cout.setf(ios::fixed);

            points.raw_print(cout, "POINTS=");
        //cout << "POINTS "<<points;
        }

        //step 2 - calculate center of gravity of all points except largest energy
        calculateCenterOfGravityDM();
        //step 3 - calculate reflection point
        calculateReflectedPointDM();
        //step 4 - decide case dependent on the energy of the reflected point
        decideCaseDM();

        //deciding to stop the algorithm
        averageEnergy=0.0;
        for(int i=0; i<points.n_rows;i++){
            currentAverageAlpha+=points(i,0);
            currentAverageBeta+=points(i,1);
            currentAverageDist+=points(i,2);
        }
        currentAverageAlpha/=points.n_rows;
        currentAverageBeta/=points.n_rows;
        currentAverageDist/=points.n_rows;

        difference=sqrt(pow((currentAverageAlpha-previousAverageAlpha),2)+pow((currentAverageBeta-previousAverageBeta),2)+pow((currentAverageDist-previousAverageDist),2));

//        if(difference<tolerance){
//            cout << "I am rank "<<myrank<<" and I reached tolerance level "<< difference<<endl;
//        }
        previousAverageAlpha=currentAverageAlpha;
        previousAverageBeta=currentAverageBeta;
        previousAverageDist=currentAverageDist;
        cout << "finishing iteration "<<endl;
    }
    cout << "out of loop "<<endl;
}

void Minimizer::calculateCenterOfGravityDM(){
    gravityalpha=0.0;
    gravitybeta=0.0;
    gravitydist=0.0;
    for(int i=0; i<points.n_rows-1; i++){
        gravityalpha+=points(i,0);
        gravitybeta+=points(i,1);
        gravitydist+=points(i,2);
    }
    gravityalpha/=(points.n_rows-1);
    gravitybeta/=(points.n_rows-1);
    gravitydist/=(points.n_rows-1);
}

void Minimizer::calculateReflectedPointDM(){
    reflectedalpha=(1.0+rhoNM)*gravityalpha-rhoNM*points(points.n_rows-1,0);
    reflectedbeta=(1.0+rhoNM)*gravitybeta-rhoNM*points(points.n_rows-1,1);
    reflecteddist=(1.0+rhoNM)*gravitydist-rhoNM*points(points.n_rows-1,2);

    if(reflectedbeta<0.0){
        reflectedbeta=-reflectedbeta;
    }
    vmcsolver.setAlphaBeta(reflectedalpha,reflectedbeta);
    vmcsolver.setdist(reflecteddist);
    energyReflectedPoint=vmcsolver.runMonteCarloIntegration();
}

void Minimizer::expansionDM(){
    expansionalpha=gravityalpha-chiNM*(gravityalpha-reflectedalpha);
    expansionbeta=gravitybeta-chiNM*(gravitybeta-reflectedbeta);
    expansiondist=gravitydist-chiNM*(gravitydist-reflecteddist);
    if(expansionbeta<0.0){
        expansionbeta=-expansionbeta;
    }
    vmcsolver.setAlphaBeta(expansionalpha,expansionbeta);
    vmcsolver.setdist(expansiondist);
    energyExpansionPoint=vmcsolver.runMonteCarloIntegration();

    if(energyExpansionPoint<energyReflectedPoint){
        if(myrank==0)
        cout <<"accepting expansion "<<endl;
        points(points.n_rows-1,0)=expansionalpha;
        points(points.n_rows-1,1)=expansionbeta;
        points(points.n_rows-1,2)=expansiondist;
        points(points.n_rows-1,3)=energyExpansionPoint;
    }
    else{
        acceptReflectionDM();
    }
}

void Minimizer::acceptReflectionDM(){
    if(myrank==0)
    cout << "accept reflection "<<endl;
    points(points.n_rows-1,0)=reflectedalpha;
    points(points.n_rows-1,1)=reflectedbeta;
    points(points.n_rows-1,2)=reflecteddist;
    points(points.n_rows-1,3)=energyReflectedPoint;
}

void Minimizer::innercontractionDM(){

    contractionalpha=gravityalpha-gammaNM*(gravityalpha-points(points.n_rows-1,0));
    contractionbeta=gravitybeta-gammaNM*(gravitybeta-points(points.n_rows-1,1));
    contractiondist=gravitydist-gammaNM*(gravitydist-points(points.n_rows-1,2));
    if(contractionbeta<0.0){
        contractionbeta=-contractionbeta;
    }
    vmcsolver.setAlphaBeta(contractionalpha,contractionbeta);
    vmcsolver.setdist(contractiondist);
    energyContractionPoint=vmcsolver.runMonteCarloIntegration();

    //accepting contraction point
    if(energyContractionPoint<points(points.n_rows-1,3)){
        if(myrank==0)
        cout << "accept innercontraction "<<endl;
        points(points.n_rows-1,0)=contractionalpha;
        points(points.n_rows-1,1)=contractionbeta;
        points(points.n_rows-1,2)=contractiondist;
        points(points.n_rows-1,3)=energyContractionPoint;
    }
    else{
        reductionDM();
    }
}

void Minimizer::outercontractionDM()
{
    contractionalpha=gravityalpha+gammaNM*(reflectedalpha-gravityalpha);
    contractionbeta=gravitybeta+gammaNM*(reflectedbeta-gravitybeta);
    contractiondist=gravitydist+gammaNM*(reflecteddist-gravitydist);
    if(contractionbeta<0.0){
        contractionbeta=-contractionbeta;
    }
    vmcsolver.setAlphaBeta(contractionalpha,contractionbeta);
    vmcsolver.setdist(contractiondist);
    energyContractionPoint=vmcsolver.runMonteCarloIntegration();

    //accepting contraction point
    if(energyContractionPoint<points(points.n_rows-1,3)){
        if(myrank==0)
        cout << "accept outercontraction "<<endl;
        points(points.n_rows-1,0)=contractionalpha;
        points(points.n_rows-1,1)=contractionbeta;
        points(points.n_rows-1,2)=contractiondist;
        points(points.n_rows-1,3)=energyContractionPoint;
    }
    else{
        reductionDM();
    }
}

void Minimizer::reductionDM(){
    if(myrank==0)
    cout << "reduction "<<endl;
    for(int i=1; i<points.n_rows;i++){
        points(i,0)=points(0,0)+sigmaNM*(points(i,0)-points(0,0));
        points(i,1)=points(0,1)+sigmaNM*(points(i,1)-points(0,1));
        points(i,2)=points(0,2)+sigmaNM*(points(i,2)-points(0,2));
        if(points(i,1)<0){
            points(i,1)=-points(i,1);
        }
        vmcsolver.setAlphaBeta(points(i,0),points(i,1));
        vmcsolver.setdist(points(i,2));
        double reductedEnergy=vmcsolver.runMonteCarloIntegration();
        points(i,points.n_cols-1)=reductedEnergy;
    }
}
