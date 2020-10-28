#ifndef CORRFILEREADER_H
#define CORRFILEREADER_H

#include "Defs.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>


namespace BiometricCorrelations{
    struct biometriccorr{
        double mean;
        double stdev;
        std::vector<double> corr_alpha;
        std::vector<double> corr_beta;
    };
    
    typedef  std::map<std::string, biometriccorr*> BiometricCorrMapType;
    
    BiometricCorrMapType ReadCorrelationFile(std::string filename); 
}

#endif
