#include "GUI/CorrelationFileReader.h"

namespace BiometricCorrelations{
    
     BiometricCorrMapType ReadCorrelationFile(std::string filename){
    
        // https://www.reddit.com/r/cpp/comments/2e68nd/stdstod_is_locale_dependant_but_the_docs_does_not/
        std::setlocale(LC_ALL, "C"); // C uses "." as decimal-point separator
         
        BiometricCorrMapType biomap;
        
        std::ifstream fin(filename.c_str());
        std::string line;
        
        if(fin.is_open()){
            cout<<"file is open"<<endl;
        } else{
            cout<<"file isn't open"<<endl;
        }

        std::string metricname = "";
        
        while(std::getline(fin, line, '\n')) {
            std::string phrase;
            std::vector<std::string> row;
            std::stringstream ss(line);
            while(std::getline(ss, phrase, '\t')) {
                row.push_back(std::move(phrase));
            }
            
            if(row.size()==0){ continue;}
            
            if(row[0]=="Metric:"){
                metricname = row[1];
                biomap[metricname] = new biometriccorr();
            }
            else if(row[0]=="Mean:"){
                biomap[metricname]->mean = std::stod(row[1]);
            }
            else if(row[0]=="Stdev:"){
                biomap[metricname]->stdev = std::stod(row[1]);
            }
            else if(row[0]=="Regression coeff alpha"){
                for(int i=1; i<row.size(); i++){biomap[metricname]->corr_alpha.push_back(std::stod(row[i]));}
            }
            else if(row[0]=="Regression coeff beta"){
                for(int i=1; i<row.size(); i++){biomap[metricname]->corr_beta.push_back(std::stod(row[i]));}
            }
        }
        return biomap;
    }    
    
    
    
}
