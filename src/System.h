#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Core>
using namespace std;
using Eigen::MatrixXd;
class System{
  private: 
    string name;
    string filename;
    int numobservables;
    MatrixXd data;
  public:
    System(string name_,string filename_,int n_);
    void setName(string name_);
    void setFileName(string fname_);
    void setNumObservables(int n_);
    string getFileName();
    string getName();
    int getNumObservables();
    void readData();
    void getData(const MatrixXd &data_);
    int getDimensionality();
};
#endif
