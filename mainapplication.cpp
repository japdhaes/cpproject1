#include "mainapplication.h"
using namespace std;

MainApplication::MainApplication(int _myrank, int _numprocs):
    myrank(_myrank),
    numprocs(_numprocs),
    alphaHe(1.819),
    betaHe(0.352),
    alphaBe(3.980),
    betaBe(0.087),
    alphaNe(10.24366653711),
    betaNe(0.09522905708),
    alphaH2(1.382),
    betaH2(0.208),
    distH2(1.263),
    alphaBe2(3.84229201074),
    betaBe2(0.17401323480),
    distBe2(5.26929149735)
{
}


void MainApplication::runApplication(){
    //==========================
    //optimal He parameters
    //energy
    //alpha = 1.819
    //beta = 0.352
    //variance
    //alpha = 1.957
    //beta = 0.361
    //==========================
    //==========================
    //optimal Be parameters
    //energy
    //alpha = 3.980
    //beta = 0.087
    //variance
    //alpha = 3.953
    //beta = 0.655
    //==========================
    //==========================
    //optimal Ne parameters
    //energy
    //alpha = 10.24366653711
    //beta = 0.09522905708
    //variance
    //alpha =10.030
    //beta = 0.350
    //==========================
    //==========================
    //optimal H2 parameters
    //energy
    //alpha = 1.382
    //beta = 0.208
    //dist = 1.263
    //==========================

    //runSimulation();
    simulateWithOutput(8, alphaBe2, betaBe2);
    //minimizeBruteForce();
    //minimizeNM();
}



void MainApplication::minimizeBruteForce(){
    Minimizer minimizer(myrank,numprocs, 4);
    double alphamin = 3.5;
    double alphamax = 4.0;
    double betamin = 0.0;
    double betamax = 0.5;
    minimizer.bruteForce(alphamin, alphamax, betamin, betamax);
}

void MainApplication::minimizeNM(){
    Minimizer minimizer(myrank, numprocs, 8);
    minimizer.nelderMeadMethodDM();
}

void MainApplication::simulateWithOutput(int nParticles, double alpha, double beta){
    TestDirectory test;
    ostringstream outputfolder;
    outputfolder << "/home/jonathan/projectsCP/project1/data/";
    //creating outputdirectory
    if(!test.DirectoryExists(outputfolder.str().c_str())){
        cout << "Output directory doesn't exist."<<endl;
        cout << "Creating output directory "<<outputfolder.str()<<endl;
        ostringstream systemmessage;
        systemmessage << "mkdir "<<outputfolder.str();
        system(systemmessage.str().c_str());
    }
    //test if creating outputdirectory worked
    if(!test.DirectoryExists(outputfolder.str().c_str())){
        cout << "Error in creating output directory. ";
        cout << "Change output directory in the script and make sure that the program ";
        cout << "is having the correct permissions to create the output directory. "<<endl;
        exit(1);
    }
    //test number of particles
    if(nParticles !=2 && nParticles != 4 && nParticles != 10){
        cout << "The program is being run with "<<nParticles<< " particles. ";
        cout << "The program has never been tested with this number of particles, ";
        cout << "so the produced data may contain errors. "<<endl;
    }
    //creating outputfolder for data
    if(nParticles == 2){
        outputfolder << "dihydrogen/";
    }
    else if(nParticles == 4){
        outputfolder << "beryllium/";
    }
    else if (nParticles == 10){
        outputfolder << "neon/";
    }
    else if (nParticles ==8){
        outputfolder << "diberyllium/";
    }
    else{
        outputfolder << nParticles<<"particles/";
    }
    if(!test.DirectoryExists(outputfolder.str().c_str())){
        cout << "Output directory for specific atom doesn't exist."<<endl;
        cout << "Creating output directory "<<outputfolder.str()<<endl;
        ostringstream systemmessage;
        systemmessage << "mkdir "<<outputfolder.str();
        system(systemmessage.str().c_str());
    }
    //test if creating outputdirectory worked
    if(!test.DirectoryExists(outputfolder.str().c_str())){
        cout << "Error in creating output directory. ";
        cout << "Change output directory in the script and make sure that the program ";
        cout << "is having the correct permissions to create the output directory. "<<endl;
        exit(1);
    }
    //Beryllium
//    int nParticles=10;
//    double alpha=10;
//    double beta=0.2;

    //Neon
//    int nParticles=10;
//    double alpha=10.6;
//    double beta=0.1;

    VMCIS vmc=VMCIS(this->myrank, this->numprocs, nParticles, alpha, beta);
    vmc.setdist(distBe2);
    vmc.setCycles(1e6);
    vmc.setOutput(true);
    vmc.setOutputDirectory(outputfolder.str().c_str());
    vmc.runMonteCarloIntegration();
}

void MainApplication::runSimulation(bool importancesampling, int nCycles, int nParticles, double alpha, double beta){
    if(importancesampling){
        VMCIS vmc=VMCIS(this->myrank, this->numprocs, nParticles, alpha, beta);
        vmc.setCycles(nCycles);
        vmc.runMonteCarloIntegration();
    }
    else{
        VMCBF vmc=VMCBF(this->myrank, this->numprocs, nParticles, alpha, beta);
        vmc.setCycles(nCycles);
        vmc.runMonteCarloIntegration();
    }
}

void MainApplication::runSimulation(){
    int nCycles=1e6;
    //Helium
    int nParticles=8;
    double alpha=3.76;
    double beta=0.49;

    //Beryllium
//    int nParticles=10;
//    double alpha=3.54406;
//    alpha=10.0;
//    double beta=0.476959;

    //Neon
//    int nParticles=10;
//    double alpha=9.8;
//    double beta=0.4;

    VMCIS vmc=VMCIS(this->myrank, this->numprocs, nParticles, alpha, beta);
    vmc.setdist(2.45);
    vmc.setCycles(nCycles);
    cout << vmc.runMonteCarloIntegration()<<endl;
}

//===============================
//FIND a solution for the fact that getAcceptanceRatio a "double (VMCBF::*)(double)" function is, while rtsec needs a
//double (*)(double) function...
//===============================
//void MainApplication::steplengthSecant(){
//    //double rtsec(double (*func)(double), double x1, double x2, double xacc)
//    VMCBF vmcbf(myrank, numprocs, 2, alphaHe, betaHe);
//    double steplength = rtsec(&vmcbf.getAcceptanceRatio,0,5,0.01);
//}

