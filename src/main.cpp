#include <string>
#include <csignal>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <LBFGS.h>
#include <fstream>
#include "lepton/Lepton.h"

#include "CustomLeptonConstants.h"
#include "System.h"
#include "Print.h"
#include "ReadFromFile.h"
#include "ObjFunction.h"
using namespace std;
//using Eigen::VectorXd;
using namespace Eigen;
using namespace LBFGSpp;

namespace
{
    volatile std::sig_atomic_t gSignalStatus;
}

void signal_handler(int signal)
{
    gSignalStatus = signal;
}


class  CombineSystems{
  private:
    vector<string> functions;
    std::vector< std::vector <string > > functions_var;
    std::vector< std::vector <unsigned > > functions_var_index;
    string lambdafilename;
    std::vector<Lepton::CompiledExpression> expression;
    std::vector< std::vector< Lepton::CompiledExpression > > expression_deriv; //Likelihood derivatives wrt observable
    vector<string> var;
    vector<System> Systems;
    vector<ObjF> SystemGradients;
    vector<double> lambda_; //if initial guess for lambda is read from input file 
    vector<ofstream> printout;
    int dimensionality;
    int maxiter;
    int leaveoutindex;
    double epsilon;
    double alpha; //penalty factor on Lagrangian multipliers
    vector<Print> Prints;
    bool leaveoneout;
    bool debug;
    bool do_minimize;
    vector<double> fit_from,fit_to;
    vector<pair<double,double>> fit_boundaries;
    std::vector< std::string > reserved_keywords = {"total_cost"};
  public:
    
    CombineSystems(){     
      //Some default parameters
      leaveoneout = false;
      debug = false;
      epsilon = 1e-6;
      maxiter = 10;
      dimensionality = -1;
      alpha=0.0;
      lambda_.resize(0);
      lambdafilename = "";
      leaveoutindex = 0; //disabled
      do_minimize =true;
      ReadFromFile(Systems,functions,maxiter,epsilon,lambda_,alpha,Prints,lambdafilename,leaveoneout,debug,do_minimize,fit_from,fit_to);
      //Construct a vector of pair with start,end of each block to use as training set
      if(fit_from.size() != fit_to.size()){
        cout<<"ERROR: Fit_from and fit_to must have the same dimension"<<endl;
        exit(1);
      }else{
        if(fit_from.size() == 0){
          fit_boundaries.push_back(make_pair(0.,100.));
        }else{
          for(int i=0;i<fit_from.size();i++)
            fit_boundaries.push_back(make_pair(fit_from[i],fit_to[i]));
        }
      }
      functions_var.resize(functions.size());
      functions_var_index.resize(functions.size());
      ofstream lambdaout;
      lambdaout.open(lambdafilename);
      cout<<"LOG: Parsed input file:"<<endl;      
      //create the variables vector. If more than one observable is present in each system, variables will be named as "name"+"0,1,2,3.."
      if(Systems.size() == 0) {
        cout<<"ERROR: No system found"<<endl;
        exit(EXIT_FAILURE);
      }
      dimensionality=Systems[0].getDimensionality();
      if(debug) cout<<"LOG: Checking dimensionality"<<endl;
      bool check_dim_failed=false;
      for(int i=0;i<Systems.size();i++){
        if(Systems[i].getDimensionality() != dimensionality)
          check_dim_failed=true;
      }
      if(check_dim_failed){
        cout<<"ERROR: Systems dimensionality are not the same"<<endl;
        cout<<"Hint: Forgot bias column?"<<endl;
        cout<<"Hint: Check the number of observable you pass to add_system command"<<endl;
        cout<<"Wrong trajectory file? N.B. The file should contain an header specifying number of rows and number of columns"<<endl;
        exit(EXIT_FAILURE);
      }
      if(debug) cout<<"LOG: Starting loop over systems"<<endl;
      for(int i=0;i<Systems.size();i++){
        cout<<"LOG: System "<<i+1<<" "<<Systems[i].getName()<<" "<<Systems[i].getFileName()<<" "<<Systems[i].getNumObservables()<<endl;
        if(debug) cout<<"LOG: Calling gradient routine constructor"<<endl;
        ObjF SystemGradient(Systems[i].getFileName(),Systems[i].getNumObservables(),fit_boundaries);
        SystemGradients.push_back(SystemGradient);
        int n=Systems[i].getNumObservables();
        if(n>1){
          for(int j=0;j<n;j++){
            string iobs=to_string(j);
            var.push_back(Systems[i].getName()+ iobs); //this should be parsed in input file
          }
        }else{
        var.push_back(Systems[i].getName()); //this should be parsed in input file
        }
      }
      cout<<"LOG: Epsilon: "<<epsilon<<endl;
      cout<<"LOG: Max. Iterations: "<<maxiter<<endl;
      cout<<"LOG: Leaveoneout: "<<leaveoneout<<endl;
      int nlike = functions.size();
      cout<<"LOG: "<<nlike<<" Likelihood functions found;"<<endl;
      
      //Parse print functions
      if(Prints.size() >0){
          cout<<"Print functions found:"<<endl;
          for(unsigned i=0;i<Prints.size();i++){
            cout<<"Print #"<<i+1<<" on file "<<Prints[i].filename<<":\n";
            std::vector<string> observables = Prints[i].observables;
            ofstream outfile;
            outfile.open(Prints[i].filename); //open the file just to cancel previously generated ones
            printout.push_back(std::move(outfile)); //store output stream
            for(unsigned j=0;j<observables.size();j++){

              if(std::find(reserved_keywords.begin(),reserved_keywords.end(),observables[j])==reserved_keywords.end()){ //check if printing a predefined variables
                Lepton::ParsedExpression pe=Lepton::Parser::parse(observables[j]).optimize(leptonConstants);
                Prints[i].LeptonCompiledExpressions.push_back(pe.createCompiledExpression());
                cout<<"\t"<<pe<<endl;
                Lepton::CompiledExpression expression = Prints[i].LeptonCompiledExpressions[j];
                for(auto &p: expression.getVariables()) {
                  if(std::find(var.begin(),var.end(),p)==var.end()) {
                    cout<<"variable " + p + " is not defined"<<endl;
                    exit(0);
                  }
                } 
              }
            }
          
        }
      }
      expression_deriv.resize(functions.size()); //as many rows as the numbers of likelihood functions
      expression.resize(functions.size());
			
      cout<<"Var size: "<<var.size()<<endl;
			//#pragma omp parallel for
      for(unsigned ifunc=0; ifunc<functions.size();ifunc++){
        //Parse functions    
        expression_deriv[ifunc].resize(var.size());
        string func = functions[ifunc];
        Lepton::ParsedExpression pe=Lepton::Parser::parse(func).optimize(leptonConstants);
        if(debug) cout<<"  Function "<<ifunc<<" as parsed by lepton: "<<pe<<"\n";
        expression[ifunc]=pe.createCompiledExpression();
        //Check if in some functions is presente a variable not defined
				for(auto &p: expression[ifunc].getVariables()) {
        if(std::find(var.begin(),var.end(),p)==var.end()) {
          cout<<"variable " + p + " is not defined"<<endl;
          exit(0);
          }
        }
        if(debug)
          cout<<"Function "<<ifunc<<"  derivatives as computed by lepton:\n";
        std::vector<string> ifuncvar;
        std::vector<unsigned> ifuncvar_index;
        for(unsigned i=0; i<var.size(); i++) {
          auto pe=Lepton::Parser::parse(func).optimize(leptonConstants);
          auto ce = pe.createCompiledExpression();
          bool found=false;
          for(auto& p : ce.getVariables()){
            if(p==var[i])found=true;
          }
          if(found){
            pe=pe.differentiate(var[i]);
            ifuncvar.push_back(var[i]);
            ifuncvar_index.push_back(i);    
          }else {
            pe=Lepton::Parser::parse("0.0");
          }
          if(debug) cout<<"w.r.t.  "<<var[i]<<" =  "<<pe<<"\n";
          expression_deriv[ifunc][i]=pe.createCompiledExpression(); //Vector of Likelihood derivatives w.r.t. each obeservable
          }
        functions_var[ifunc]=ifuncvar;
        functions_var_index[ifunc]=ifuncvar_index;
        if(debug) cout<<"LOG: Derivatives computed correctly"<<endl;
        
      }
    }; //end constructor
    int getMaxIter(){
      return maxiter;
    }
    int getDimensionality(){
      return dimensionality;
    }
    double getEpsilon(){
      return epsilon;
    }
    int getNumLikelihoods(){
      return functions.size();
    }
    void setLeaveOut(int indexout){
      leaveoutindex = indexout;
    }
    bool getLeaveOneOut(){
      return leaveoneout;
    }
    bool getDoMinimize(){
      return do_minimize;
    }

    bool firstiter;
    double totlikelihood;
    double totlikelihood_0;
    void CheckLambdaFromFile(VectorXd& lambda){
      if(lambda_.size()>0){
        if(lambda_.size()!=dimensionality){
          cout<<"ERROR: Lambda dimensions from input file don't match dimensionailty;"<<endl;
          cout<<"The dimension of parsed vector is "<<lambda_.size()<<endl;
          cout<<"Parsed values are:"<<endl;
          for(int i=0;i<lambda_.size();i++)
            cout<<lambda_[i]<<endl;
          exit(EXIT_FAILURE);
        }else{
          for(int i=0;i<lambda.rows();i++)
            lambda[i]=lambda_[i];
          cout<<"LOG: Reading Lagrangian Multipliers from file; Values are:\n";
          cout<<lambda.transpose()<<endl;
        }
      }
    }
    double operator()(const VectorXd& lambda, VectorXd& grad_likelihood){
     
		  //cout<<"Start operator"<<endl;
      std::signal(SIGINT, signal_handler);
      if(gSignalStatus == 2){
        cout<<"Stopping ... "<<endl;
        cout<<"lambdas = \n"<<lambda.transpose()<<endl;
        exit(1);
      }

        
      MatrixXd observables_deriv(0,dimensionality);

			VectorXd observables;
      for(int i=0;i<Systems.size();i++){

        VectorXd observables_;
        MatrixXd observables_deriv_;
        VectorXd system_obs; //observables of the ith system
				MatrixXd system_obs_deriv; //derivatives of the ith system

        observables_=observables; //copy observables computed up to now in a tmp vec
        observables_deriv_=observables_deriv; //same as before
        

				//cout<<"Calling gradient routine"<<endl;
        if(debug) cout<<"LOG: Calling gradient routine"<<endl;
        SystemGradients[i](system_obs,system_obs_deriv,lambda);
        
        observables.resize(observables.rows()+system_obs.rows()); //resize observables to host new observables
        observables_deriv.resize(observables_deriv.rows()+system_obs_deriv.rows(),dimensionality); //same as before
				
				observables<<observables_,system_obs; //concatenate to observables the ones computed on each system
				observables_deriv<<observables_deriv_,system_obs_deriv;
        
			}
		
			//cout<<"Observables deriv: (should be Nobs x dimensionality) "<<observables_deriv<<endl;

			//Evaluate Likelihood derivatives
			
      double funcval = 0.0;
      totlikelihood= 0.0;
		  grad_likelihood.fill(0.0);
     
      VectorXd likelihood_deriv; //Evaluated Likelihood derivatives wrt to observables
      likelihood_deriv.resize(var.size());
      for(unsigned ifunc=0;ifunc<functions.size();ifunc++){
        likelihood_deriv.fill(0.0);
        for(int i=0;i<functions_var_index[ifunc].size();i++){
          int ivar=functions_var_index[ifunc][i];
          try{expression[ifunc].getVariableReference(var[ivar])=observables[ivar];}catch(Lepton::Exception& exc){}
          for(int j=0;j<functions_var_index[ifunc].size();j++)
            try{
              int jvar=functions_var_index[ifunc][j];
              expression_deriv[ifunc][ivar].getVariableReference(var[jvar])=observables[jvar];
            }catch(Lepton::Exception& exc){}
          likelihood_deriv[ivar]=expression_deriv[ifunc][ivar].evaluate();
        }
        if(ifunc!=leaveoutindex){
          grad_likelihood += 1.0*(( ( likelihood_deriv.transpose() ) * observables_deriv ).transpose());
          funcval+=expression[ifunc].evaluate();
        }
        totlikelihood+=expression[ifunc].evaluate(); //accumulate total likelihood
      }
      double penalty=0.0;
      penalty=0.5*alpha*lambda.dot(lambda);
      grad_likelihood+=alpha*lambda;

      /*
      for(int i=0;i<lambda.rows();i++){
        penalty+=lambda[i]*lambda[i];
        grad_likelihood[i]+=alpha*lambda[i];
      }
      penalty*=0.5*alpha;
*/

      for(unsigned i=0;i<Prints.size();i++){
        ofstream out;
        out.open(Prints[i].filename,ios::app);
          //ofstream out=std::move(printout[i]);
          for(unsigned j=0; j<var.size();j++){
            for(unsigned k=0;k<Prints[i].LeptonCompiledExpressions.size();k++){
              try{Prints[i].LeptonCompiledExpressions[k].getVariableReference(var[j])=observables[j];}catch(Lepton::Exception& exc){}
            }
          }
          for(unsigned j=0; j<Prints[i].LeptonCompiledExpressions.size();j++){
            out<<Prints[i].LeptonCompiledExpressions[j].evaluate()<<" ";
          }
        
          for(int printobs = 0;printobs < Prints[i].observables.size();printobs++){
            if(Prints[i].observables[printobs].compare("total_cost")==0)
              out<<totlikelihood+penalty<<" ";
          }
        
          out<<endl;
      }
      if(lambdafilename!=""){
        ofstream out;
        out.open(lambdafilename,ios::app);
        out<<lambda.transpose()<<endl;
      }
        
      if(debug) cout<<"Function value: "<<funcval<<"; with penalty = "<<funcval+penalty<<endl;
      if(firstiter) {totlikelihood_0=funcval; firstiter=false;}

			return funcval+penalty;
		  //cout<<"End operator"<<endl;
		}

    ~CombineSystems(){
      for(unsigned i=0;i<printout.size();i++)
        printout[i].close();
    }

};


int main()
{
  
 
 // Set up parameters
 LBFGSParam<double> param;
 CombineSystems systems;
 const int n = systems.getDimensionality(); 
 param.epsilon = systems.getEpsilon();
 param.max_iterations = systems.getMaxIter();
 int numlikelihoods = systems.getNumLikelihoods();
 int niterations;
 bool leaveoneout = systems.getLeaveOneOut();
 bool do_minimize = systems.getDoMinimize();
 if(!leaveoneout)  
   niterations = 1;
 else
   niterations = numlikelihoods;
 // Create solver and function object
 LBFGSSolver<double> solver(param);
 //Rosenbrock fun(n);
 
 
// objf myf;
 // Initial guess
 VectorXd lambda = VectorXd::Zero(n);
 VectorXd fake_grad = VectorXd::Zero(n);
 //Check if Lambda starting values are read from the input file


 for(int outindex = 0;outindex <= niterations;outindex++){
  
  systems.setLeaveOut(outindex-1);
  lambda.fill(0.0);
  systems.CheckLambdaFromFile(lambda);
  if(!do_minimize){
    systems(lambda,fake_grad);
    return 0;
  }

 // x will be overwritten to be the best point found

//VectorXd g(n);
//  lambda[0]=5;
//  std::cout<<myf(lambda,g)<<"\n";
//	return 0;

//std::cout<<systems(lambda,g)<<"\n";
 std::vector<Lepton::CompiledExpression> likelihoods;
 double funcval;
 systems.firstiter=true;
 int niter = solver.minimize(systems, lambda, funcval);
 //std::cout << niter << " iterations" << std::endl;
 //std::cout << "lambda = \n" << lambda.transpose() << std::endl;
 //std::cout << "Likelihood= " << funcval << std::endl;
 std::cout<<"\t"<<niter<<" iter.; Initial Likelihood = "<<systems.totlikelihood_0<<"; Total likelihooud leaving out "<<outindex-1<<" = "<<systems.totlikelihood<<endl;
  return 0;
}
}
