#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <algorithm>



int main(){

std::ifstream breach_linelist("scenario1/arrivals_linelist_scenario1_sim1.csv");

std::vector<double> total_time;
std::vector<double> FOI; 
std::vector<double> exposure_days;

if(breach_linelist.is_open()){
        std::string line;
        double value;
        std::getline(breach_linelist,line); // Get first line. 
        // Do nothing this is a title.
        // Secondline onwards. 
        while(std::getline(breach_linelist,line)){
            std::stringstream stream_line(line);
            std::string row_val;
            // std::vector<double> row;
            for(int col_number = 0; col_number < 6; col_number++){
                std::getline(stream_line,row_val,','); // Row val is the corresponding double we are after. 
                std::stringstream stream_row(row_val);
                stream_row >> value;
                if(col_number ==1){
                    exposure_days.push_back(value);
                }
                if(col_number ==2){
                    FOI.push_back(value);
                }
                if(col_number ==5){
                    total_time.push_back(value);
                }
            }
        }
        breach_linelist.close();
    } else {
        throw std::logic_error("The schedule file for the scenario \n");
    }

    auto printstuff = [](std::vector<double> vec)->void{
        for(auto x: vec){
            std::cout << x << ", ";
        }
        std::cout << std::endl;
    };

    printstuff(total_time);
    printstuff(FOI);
    printstuff(exposure_days);

    return 0;
}