#ifndef PRINT_H
#define PRINT_H

#include <iostream>
#include <string>
#include "lepton/Lepton.h"

using namespace std;
class Print{
  public:
    std::string filename;
    std::vector<string> observables;
    std::vector<Lepton::CompiledExpression> LeptonCompiledExpressions;
};

#endif
