#ifndef OBJFUNCTION_H
#define OBJFUNCTION_H
#include <string>
#include <iostream>
//#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using namespace Eigen;
//using Eigen::VectorXd;
//using Eigen::MatrixXd;

class ObjF{
  private:
    int dimensionality;
    int numobservables;
    int nframes;
    string datafilename;
    MatrixXd datamatrix,observables;
    VectorXd bias;
    void readData(MatrixXd &datamatrix,MatrixXd &observables,VectorXd &bias,int &nframes,const vector<pair<double,double> > &fit_boundaries);
  public:
    ObjF(string datafilename_ ,int n_,const vector<pair<double,double> > &fit_boundaries);
    void operator()(VectorXd& obs, MatrixXd& deriv,const VectorXd& lambda);
    void gradient(VectorXd& obs, MatrixXd& deriv,const VectorXd& lambda);
  //Funzione che dato il nome dell'osservabile ritorna il gradiente di quell'osservabile
  //~ObjF();
};
#endif
