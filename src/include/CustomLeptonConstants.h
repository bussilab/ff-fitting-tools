#ifndef CUSTOM_LEPTON_CONSTANTS_H_
#define CUSTOM_LEPTON_CONSTANTS_H_

#define _USE_MATH_DEFINES

#include <string>
#include <cmath>
#define pi M_PI

//using namespace std;
static std::map<std::string, double> leptonConstants= {
  {"e", std::exp(1.0)},
  {"log2e", 1.0/std::log(2.0)},
  {"log10e", 1.0/std::log(10.0)},
  {"ln2", std::log(2.0)},
  {"ln10", std::log(10.0)},
  {"pi", pi},
  {"pi_2", pi*0.5},
  {"pi_4", pi*0.25},
//  {"1_pi", 1.0/pi},
//  {"2_pi", 2.0/pi},
//  {"2_sqrtpi", 2.0/std::sqrt(pi)},
  {"sqrt2", std::sqrt(2.0)},
  {"sqrt1_2", std::sqrt(0.5)}
};



#endif /*CUSTOM_LEPTON_CONSTANTS_H_*/
