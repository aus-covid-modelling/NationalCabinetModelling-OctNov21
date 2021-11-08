#include "Quarantine_outputs.h"
#include <iostream>
#include <fstream>
#include <sstream>

void load_linelist(std::vector<double> & total_time, std::vector<double> & exposure_days, std::vector<double>& FOI, std::vector<int> & kid,std::string filename){
    
    std::ifstream breach_linelist(filename);
    if(breach_linelist.is_open()){
        std::string line;
        double value;
        std::getline(breach_linelist,line); // Get first line. 
        // Do nothing this is a title.
        // To future proof we should check if strings match what we are after, otherwise throw an error? 
        // Secondline onwards. 
        while(std::getline(breach_linelist,line)){
            std::stringstream stream_line(line);
            std::string row_val;
            // std::vector<double> row;
            for(int col_number = 0; col_number < 9; col_number++){
                std::getline(stream_line,row_val,','); // Row val is the corresponding double we are after. 
                std::stringstream stream_row(row_val);
                stream_row >> value;
                if(col_number ==1){
                    exposure_days.push_back(value);
                }
                if(col_number == 6){
                    FOI.push_back(value);
                }
                if(col_number ==5){
                    total_time.push_back(value);
                }
                if(col_number ==8){
                    kid.push_back(value);
                }
            }
        }
        breach_linelist.close();
    } else {
        throw std::logic_error("The schedule file for scenario " + filename + " was not found.");
    }
}


pseudo_individual::pseudo_individual(double t_start, double t_E, double FOI, int age_bracket_in, std::vector<double> & age_stratified_contacts,std::vector<double> & xi_k,int cluster_ref):start_exposure(t_start),finish_exposure(t_start + t_E),age_bracket(age_bracket_in),cluster(cluster_ref){
        // Calculate alpha here. 
        if(xi_k.size()!=age_stratified_contacts.size()){
            throw std::logic_error("xi_k does not match the contact matrix size in func: pseudo_individual::pseudo_individual.");
        }
        double susceptible_contacts =  0.0;
        for(int i =0; i < xi_k.size(); i++){
            susceptible_contacts += xi_k[i]*age_stratified_contacts[i];
        }

        alpha = FOI/(t_E * susceptible_contacts);
}