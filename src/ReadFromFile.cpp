#include "System.h"
#include "Print.h"
#include "ReadFromFile.h"
//Default parameters

void ReadFromFile(vector<System> &systems, vector<string> &functions,int &maxiter,double &epsilon, vector<double> &lambda,double &alpha, vector<Print>& Prints,string& lambdaout, bool &leaveoneout,bool &debug,bool &do_minimize, vector<double> &fit_from, vector<double> &fit_to){

  int nline=0;
  bool iscommand;
  string line,dummy;
  while(getline(cin,line)){
    if(line.c_str()[0]=='#'){
      nline++;
      continue;
    }
    else{
      string parameter="",command;
      istringstream iss(line);
      iscommand=true;
      do{
        if(iscommand){
          iss>>command;
          iscommand=false;
        }
        else{
          //iss>>parameter;
          if(command=="add_system"){
            iss>>parameter;
            string name_=parameter;
            string filename_;
            int n_;
            if(parameter!=""){
              if(!iss.eof()){
                iss>>parameter;
                filename_=parameter;
              }else{
                cout<<"Please specify filename and number of observables for system "+ name_ <<endl;
                exit(EXIT_FAILURE);
              }
              if(!iss.eof()){
                iss>>parameter;
                n_=atoi(parameter.c_str());
                System system_(name_,filename_,n_);
                systems.push_back(system_);
                iss>>dummy;
              }else{
                cout<<"Please specify number of observables for system "+name_ <<endl;
                exit(EXIT_FAILURE);
              }
            }else{
              cout<<"Please specify name filename and number of observables for system in the command add_system" <<endl;
              exit(EXIT_FAILURE);
            }
          nline++;
          continue;
          }
          if(command=="function"){
            iss>>parameter;
            functions.push_back(parameter);
            iss>>dummy;
            nline++;
            continue;
          }
          if(command=="print"){
            Print printcmd;
            iss>>parameter;
            printcmd.filename=parameter; //filename
            iss>>parameter; //string with observables separated by a comma
            parsePrintVariables(parameter,printcmd);
            Prints.push_back(printcmd);
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="maxiter"){
            iss>>parameter;
            maxiter=atoi(parameter.c_str());
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="epsilon"){
            iss>>parameter;
            epsilon=atof(parameter.c_str());
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="fit_from"){
            iss>>parameter;
            SplitStringToVectorF(parameter,fit_from);
            if(fit_from.size()==0){
              cout<<"ERROR: Missing parameters for fit_from."<<endl;
              exit(1);
            }
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="fit_to"){
            iss>>parameter;
            SplitStringToVectorF(parameter,fit_to);
            if(fit_to.size()==0){
              cout<<"ERROR: Missing parameters for fit_to."<<endl;
              exit(1);
            }
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="lambda"){
            while(iss>>parameter)
              lambda.push_back(atof(parameter.c_str()));
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="alpha"){
            iss>>parameter;
            alpha=atof(parameter.c_str());
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="lambdafile"){
            iss>>parameter;
            lambdaout=parameter;
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="leaveoneout"){
            leaveoneout = true;
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="debug"){
            debug = true;
            iss>>dummy;
            nline++;
          continue;
          }
          if(command=="no_minimization"){
            do_minimize = false;
            iss>>dummy;
            nline++;
          continue;
          }
        cout<<"ERROR: Unknown command "<<command<<endl;
        exit(EXIT_FAILURE);
        }
      }while(iss);
    }

  }
}

void parsePrintVariables(const std::string &str,Print& printobj){
  std::istringstream ss(str);
  std::string var;
  while(std::getline(ss, var, ';')){
    printobj.observables.push_back(var);
  }
}
void SplitStringToVectorF(const std::string &str,vector<double>& vec){
  std::istringstream ss(str);
  std::string var;
  while(std::getline(ss, var, ',')){
    vec.push_back(atof(var.c_str()));
  }
}
