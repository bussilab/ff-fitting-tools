#include "ObjFunction.h"
  ObjF::ObjF(string datafilename_ ,int n_, const vector<pair<double,double> > &fit_boundaries){
      datafilename=datafilename_;
      numobservables=n_;
			readData(datamatrix,observables,bias,nframes,fit_boundaries);
   };
  void ObjF::readData(MatrixXd &datamatrix,MatrixXd &observables,VectorXd &bias,int &nframes, const vector< pair<double,double> > &fit_boundaries){
    //fit_from and fit_to: percentage of file to use in minimization
    ifstream datafile;
    datafile.open(datafilename.c_str(),ios::binary | ios::in);
    double data;
    int ncol;
    datafile.read(reinterpret_cast<char *>(&data),sizeof(double));
    nframes=int (data);
    datafile.read(reinterpret_cast<char *>(&data),sizeof(double));
    ncol=int (data);
    dimensionality=ncol-2-numobservables;

    int ndata = 0;
    //Make a small loop to get the numbr of frames over which we are fitting
    for(int ipair=0;ipair<fit_boundaries.size();ipair++){
      double fit_from = fit_boundaries[ipair].first;
      double fit_to = fit_boundaries[ipair].second;
      int frame_start = floor(fit_from*nframes/100.);
      int frame_end =  floor(fit_to*nframes/100.);
      if(fit_from < 0 ) {
       cout<<"ERROR: fit_from can't be < 0 "<<endl;
       exit(1);
      }
      if(fit_to > 100 ) {
       cout<<"ERROR: fit_to can't be > 100 "<<endl;
       exit(1);
      }
      if(frame_end <=frame_start){
       cout<<"ERROR: frame end cannot be smaller than fram start"<<endl;
       exit(1);
      }
      ndata+=frame_end - frame_start;
    }
   
    cout<<"ndata = "<<ndata<<endl;
    datamatrix.resize(ndata,dimensionality);
    observables.resize(ndata,numobservables);
    bias.resize(ndata);
    
    MatrixXd datamatrix_all(nframes,dimensionality);
    MatrixXd observables_all(nframes,numobservables);
    VectorXd bias_all(nframes);

    //Read all frames in memory
    unsigned int i=0;
    cout<<"Reading file "<<datafilename.c_str()<<"; size "<<nframes<<endl;
    for(int iframe = 0;iframe<nframes;iframe++){
        //discard time
        datafile.read(reinterpret_cast<char *>(&data),sizeof(double));
        //read features
        for(int iarg=0;iarg<dimensionality;iarg++)
          datafile.read(reinterpret_cast<char *>(&datamatrix_all(i,iarg)),sizeof(double));
        //read observables
        for(int iarg=0;iarg<numobservables;iarg++){
          datafile.read(reinterpret_cast<char *>(&observables_all(i,iarg)),sizeof(double));
        }
        //read bias
        datafile.read(reinterpret_cast<char *>(&bias_all(i)),sizeof(double));
        i++;
    }
    datafile.close();
  

    //Read the blocks
    i = 0;
    for(int ipair=0;ipair<fit_boundaries.size();ipair++){
      double fit_from = fit_boundaries[ipair].first;
      double fit_to = fit_boundaries[ipair].second;
      int frame_start = floor(fit_from*nframes/100.);
      int frame_end =  floor(fit_to*nframes/100.);
      cout<<"Reading file "<<datafilename.c_str()<<" from frame "<<frame_start<<" to frame "<<frame_end<<endl;
      for(int iframe = frame_start;iframe<frame_end;iframe++){
        for(int iarg=0;iarg<dimensionality;iarg++){
          //cout<<"datamatrix("<<i<<","<<iarg<<") = datamatrix_all("<<iframe<<","<<iarg<<")"<<endl;
          datamatrix(i,iarg) = datamatrix_all(iframe,iarg);}
        //read observables
        for(int iarg=0;iarg<numobservables;iarg++)
          observables(i,iarg) = observables_all(iframe,iarg);
        bias(i) = bias_all(iframe);
        i++;
      }
    }
  }

   void ObjF::operator()(VectorXd& obs, MatrixXd& deriv,const VectorXd& lambda){
		VectorXd logweights = -datamatrix*lambda;
    logweights+=bias;
    double maxw = logweights.maxCoeff();
    logweights.array()-=maxw;
    VectorXd weights=(logweights.array()).exp();
    double norm=weights.sum();
    double stateff=(norm*norm)/weights.dot(weights);
    //cout<<"Stat. Eff. "<<datafilename<<": "<<100*stateff/nframes<<" %"<<endl;
    cout<<"Stat. Eff. "<<datafilename<<": "<<stateff<<endl;
    MatrixXd v1(numobservables,dimensionality);
    VectorXd rew_data(dimensionality); 
    //1. multiply each column of datamatrix for its weight 
    //2. compute Observable.transpose() * the result of 1. to get v1
    v1  = observables.transpose()*(datamatrix.array().colwise() * weights.array()).matrix();
    v1/=norm; //normalize
    obs = (observables.transpose()*weights)/norm;
    rew_data = (datamatrix.transpose()*weights)/norm;
    deriv = -v1+obs*rew_data.transpose();
   }
   
