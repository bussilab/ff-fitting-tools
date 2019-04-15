#ifndef READFROMFILE_H_INCLUDED
#define READFROMFILE_H_INCLUDED

#include <string>
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
using namespace std;

void ReadFromFile(vector<System> &systems, vector<string> &functions,int &maxiter,double &epsilon, vector<double> &lambda,double &alpha, vector<Print>& Prints,string& lambdaout,bool &leaveoneout,bool &debug, bool &do_minimize,double &fit_from, double &fit_to);
void ReadFromFile(vector<System> &systems, vector<string> &functions,int &maxiter,double &epsilon, vector<double> &lambda,double &alpha, vector<Print>& Prints,string& lambdaout, bool &leaveoneout,bool &debug,bool &do_minimize, vector<double> &fit_from, vector<double> &fit_to);
void parsePrintVariables(const std::string &str,Print& printobj);
void SplitStringToVectorF(const std::string &str,vector<double>& vec);
#endif // READFROMFILE_H_INCLUDED
