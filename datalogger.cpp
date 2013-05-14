#include "datalogger.h"

Datalogger::Datalogger(const string &outputfile, const int &lengthVector):
    i(0),
    outputFilename(outputfile),
    data(zeros(lengthVector))
{
}

void Datalogger::addData(double &energy){
    data(i)=energy;
    i++;
}

void Datalogger::flushData(){
    for(int j=0; j<data.n_rows; j++){
        (*myoutputFile) << data(j) << endl;
    }
    i=0;
}

void Datalogger::flushFinalData(){
    for(int j=0; j<i;j++){
//        myFile.write (data(j), 100);
        (*myoutputFile) << data(j) << endl;
    }
    (*myoutputFile).close();
}

void Datalogger::initialize(){
    this->myoutputFile=new ofstream();
    this->myoutputFile->open(outputFilename.c_str(), ios::out |std::ios::binary);
}

//ostringstream temp;
//temp << outputDirectory << "data"<< myrank<<".dat";
//this->outputFilename=temp.str();
//cout<<outputFilename.c_str()<<endl;

//this->myoutputFile->open(outputFilename.c_str(),std::fstream::out);
