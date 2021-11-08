#ifndef QUARANTINE_OUTPUTS_H
#define QUARANTINE_OUTPUTS_H
#include <vector>
#include <string>
void load_linelist(std::vector<double> & total_time, std::vector<double> & exposure_days, std::vector<double>& FOI, std::vector<int> & kid, std::string filename);

class pseudo_individual{
    public:
    pseudo_individual(double t_start, double t_E, double FOI, int age_bracket_in, std::vector<double> & age_stratified_contacts,std::vector<double> & xi_k,int cluster_ref);
    double start_exposure;
    double finish_exposure;
    double alpha; // Scale to ensure expected secondary infections is FOI from quarantine. 
    int age_bracket; 
    int cluster; 
};
#endif