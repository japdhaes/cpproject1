#include "mainapplication.h"

MainApplication::MainApplication(int _myrank, int _numprocs):
    myrank(_myrank),
    numprocs(_numprocs)
{
}

using namespace std;

void MainApplication::runApplication(){
    //runSimulation();
    //simulateWithOutput(10, 9.95, 0.18);
//    minimizeBruteForce();
//    minimizeNM();
    //rtsec(double (*func)(double), double x1, double x2, double xacc)

    cout << "out of minimizer "<<endl;
}



void MainApplication::minimizeBruteForce(){
    Minimizer minimizer(myrank,numprocs, 2);
    minimizer.bruteForce(10, 9, 11, 0, 1);
}

void MainApplication::minimizeNM(){
    Minimizer minimizer(myrank, numprocs, 2);
    minimizer.nelderMeadMethod();
}

void MainApplication::blockData(){
    Blocking blocking("/home/jonathan/projectsCP/project1/data/helium/");
    blocking.doBlocking();
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
        outputfolder << "helium/";
    }
    else if(nParticles == 4){
        outputfolder << "beryllium/";
    }
    else if (nParticles == 10){
        outputfolder << "neon/";
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
    bool importancesampling=true;
    //Helium
//    int nParticles=2;
//    double alpha=1.8;
//    double beta=0.36;

    //Beryllium
//    int nParticles=4;
//    double alpha=3.54406;
//    double beta=0.476959;

    //Neon
    int nParticles=10;
    double alpha=9.8;
    double beta=0.4;

    VMCIS vmc=VMCIS(this->myrank, this->numprocs, nParticles, alpha, beta);
    vmc.setCycles(nCycles);
    cout << vmc.runMonteCarloIntegration()<<endl;
}


//void MainApplication::runBeryllium(double alpha, double beta){
////    double alpha =3.5;
////    double beta=3.5;
//    BeWF solver(myrank, numprocs, alpha, beta);
//    solver.importanceSampling=true;
//    solver.closedForm=false;
//    cout << solver.runMonteCarloIntegration()<<endl;

//}

//void MainApplication::steplength_secant(){
//    double delta;

//    double steplength_pp=1, steplength_p=1.1;
//    double fpp, fp=-1, f;

//    double steplength=1;
//    double alpha=1.8;
//    double beta=0.36;

//    int varsteps=100;
//    double steplengthmin=1.4;
//    double steplengthmax=1.6;
//    double dr=(steplengthmax-steplengthmin)/varsteps;
//    double energy;
//    double tolerance=1e-4;

//    HeWF solver(myrank, numprocs, alpha, beta);
////    VMCSolver solver(alpha, beta);
//    double olderdelta, olddelta;
//    double oldersteplength=1.1, oldsteplength=1.1;

////    solver.runMonteCarloIntegration(&energy, steplength_pp, fpp);
//    fpp-=0.5;
//    //solver.runMonteCarloIntegration(&energy, steplength, delta);

//    //cout << abs(delta-0.5)/delta << endl;



//    while(abs(fpp)>tolerance){

////        solver.runMonteCarloIntegration(&energy, steplength_p, fp);
//        fp-=0.5;

//        steplength = steplength_p - fp * (steplength_p-steplength_pp)/(fp - fpp);
//        steplength_pp = steplength_p;
//        steplength_p = steplength;
//        fpp=fp;
//        cout << "steplength " << steplength << endl;
//    }

//}

//void MainApplication::steplength_manually(){
//    double steplength;
//    double alpha=1.8;
//    double beta=0.36;
//    double delta;

//    int varsteps=100;
//    double steplengthmin=1.4;
//    double steplengthmax=1.6;
//    double dr=(steplengthmax-steplengthmin)/varsteps;


//    double lowestenergy=0, energy;
//    double lowestalpha=0, lowestbeta=0;

//    HeWF solver(myrank, numprocs, alpha, beta);

////    VMCSolver solver(alpha, beta);
//    for(double nsteplength=0; nsteplength<varsteps; nsteplength++){
//        steplength = steplengthmin + nsteplength*dr;
//        if(solver.myrank==0) cout << "steplength   = " << steplength << endl;
////        solver.runMonteCarloIntegration(&energy, steplength, delta );
//    }
//}

//void MainApplication::alphabetavalues(){
//    double lowestenergy=0, energy;
//    double lowestalpha=0, lowestbeta=0;
//    /*long seed=-1*time(0);
//    cout << "test" << endl;
//    cout << "Hello World!" << endl;
//    for(int i=0; i<50; i++){
//        cout << ran2(&seed) << endl;
//    }*/


//    /*for(double i=1.78; i<1.81; i+=0.005){
//        for(double j=0.25; j<0.5; j+=0.01){

//            VMCSolver test(1, 1);
//            cout << "testing with alpha= " << i << " and beta= " <<j << endl;
//            test.runMonteCarloIntegration(&energy);
//            if(energy<lowestenergy){
//                lowestenergy=energy;
//                lowestalpha=i;
//                lowestbeta=j;
//            }
//            //found lowest energy -2.89226 at values alpha = 1.85 and beta = 0.45
//            //found lowest energy -2.89395 at values alpha = 1.8 and beta = 0.36

//        //}
//    }*/
//    cout << "found lowest energy "<< lowestenergy << " at values alpha = " << lowestalpha << " and beta = " << lowestbeta << endl;

//}

//void MainApplication::calculateClosedForm(){
//    double steplength=1.485;
//    double energy;
//    double delta;
//    double alpha=1.8;
//    double beta=0.36;

//    HeWF solver(myrank, numprocs, alpha, beta);

////    VMCSolver solver(alpha, beta);



////    solver.runMonteCarloIntegrationClosedForm(&energy, steplength, delta);

//    cout << energy << endl;
////    solver.runMonteCarloIntegration(&energy, steplength, delta);
////    cout << energy << endl;
//    //solver.runMonteCarloIntegration(&energy, steplength, delta);

//    //cout << abs(delta-0.5)/delta << endl;

//}

//void MainApplication::calculateCFImportanceSampling(){
//    double steplength=1.485;
//    double energy;
//    double delta;
//    double alpha=1.8;
//    double beta=0.36;

//    HeWF solver(myrank, numprocs, alpha, beta);

////    VMCSolver solver(alpha, beta);



////    solver.runMonteCarloIntegrationCFImportanceSampling(&energy, steplength, delta);

//    cout << energy << endl;
////    solver.runMonteCarloIntegration(&energy, steplength, delta);
////    cout << energy << endl;
//    //solver.runMonteCarloIntegration(&energy, steplength, delta);

//    //cout << abs(delta-0.5)/delta << endl;

//}
