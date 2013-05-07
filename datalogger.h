#ifndef DATALOGGER_H
#define DATALOGGER_H
#include <iostream>
#include <fstream>
#include <armadillo>
using namespace std;
using namespace arma;

class Datalogger
{
public:
    Datalogger(){};
    Datalogger(const string &outputfile, const int &lengthVector);
    void initialize();
    void addData(double &energy);
    void flushData();
    void flushFinalData();
protected:
    int             i;
    string          outputFilename;
    ofstream        *myoutputFile;
    vec             data;


};

#endif // DATALOGGER_H
