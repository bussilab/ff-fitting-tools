#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
using namespace std;

void check_bin_file(const string &name){
  double data;
  ifstream infile;
  infile.open(name.c_str(),ios::binary | ios::in);
  while(infile.read(reinterpret_cast<char *>(&data),sizeof(double))){
    //infile.read(reinterpret_cast<char *>(&data),sizeof(double));
    cout<<data<<" "<<endl;
    }
  infile.close();
}
int main(int argc, char *argv[]){
  if(argc!=3){
    cout<<"The program expect 3 arguments as follow:"<<endl;
    cout<<"file2bin.x in.txt out.bin"<<endl;
    exit(EXIT_FAILURE);
  }
  string inname=argv[1];
  string outname=argv[2];
  ifstream infile;
  infile.open(inname.c_str());
  ofstream outfile;
  outfile.open(outname.c_str(),ios::binary | ios::out);
  double data;
  if(infile.good()){
    string str;
    while(getline(infile,str)){
      istringstream ss(str);
      while(ss>>data){
        //cout<<data<<" ";
        outfile.write(reinterpret_cast<char *>(&data),sizeof(double));
      }
      //cout<<endl;
    }
  }
  outfile.close();
  infile.close();
  //check_bin_file(outname);
}

