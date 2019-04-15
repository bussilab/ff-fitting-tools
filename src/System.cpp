#include "System.h"
  
    System::System(string name_,string filename_,int n_){
      name=name_;
      filename=filename_;
      numobservables=n_;
    };
    void System::setName(string name_) {
      name=name_;
    }
    void System::setFileName(string fname_){
      filename=fname_;
    }
    void System::setNumObservables(int n_){
      numobservables=n_;
    }
      string System::getFileName(){
      return filename;
    }
      string System::getName(){
      return name;
    }
    int System::getNumObservables(){
      return numobservables;
    }
    int System::getDimensionality(){
     ifstream datafile;
     datafile.open(filename.c_str(),ios::binary | ios::in);
     double data;
     int ncol;
     datafile.read(reinterpret_cast<char *>(&data),sizeof(double));
     int nframes=int (data);
     datafile.read(reinterpret_cast<char *>(&data),sizeof(double));
     ncol=int (data);
     int dimensionality=ncol-2-numobservables;
     cout<<"System "<<name<<": "<<nframes<<"; dimensionality "<<dimensionality<<endl;
     datafile.close();
     return dimensionality;
    }

