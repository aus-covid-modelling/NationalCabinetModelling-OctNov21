#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <algorithm>
#include "../nlohmann/json.hpp"
#include "Quarantine_outputs.h"
#include "Quarantine_breach_model.h"

int main(int argc, char *argv[]){

    if(argc == 1){
        std::cout << "ERROR: Parameter file not loaded!" << std::endl;
        return 1;
    }else if(argc == 2){
        std::cout << "ERROR: Sim number not provided." << std::endl;
        return 1;
    }else if(argc ==3){
        std::cout << "ERROR: json not provided." << std::endl;
        return 1;
    }
    // Parameter read in and simulation here!
    std::string scenario_name(argv[1]); // Reads in scenario name.
    std::string json_name(argv[3]); // Reads in scenario name.
    int sim_number;

    if(argc > 4){
        std::cout << "More than 2 arguments provided to main function, ignoring extra arguments. " << std::endl;
    }
    
    // Read sim number
    std::istringstream iss(argv[2]);
    if(iss >> sim_number)
    {
        std::cout << "Sim number " << sim_number << std::endl;
    }else{
        throw std::logic_error("Sim number failed to convert\n");
    }

    // Read in scenario json file. 
    std::ifstream parameters_in(json_name);
    if(!parameters_in.is_open()){
        throw std::logic_error("Parameter json is not found. \n FILENAME : " + scenario_name + ".json\n");
    }
    // Convert input to json.
    nlohmann::json parameters_json;
    parameters_in >> parameters_json;
    parameters_in.close(); // Close the file.

    // Create the folder for outputs of the scenario. 
    std::string folder = scenario_name + (std::string) parameters_json["folder_suffix"];
    std::string directory = (std::string) parameters_json["output_directory"] + folder;

    #ifdef _WIN32
        int top_folder = mkdir(((std::string) parameters_json["output_directory"]).c_str());
        int main_folder = mkdir(directory.c_str()); // Create folder.
    #else
        int top_folder = mkdir(((std::string) parameters_json["output_directory"]).c_str(), 0777);
        int main_folder = mkdir(directory.c_str(),0777); // Create folder.
    #endif
    (void) main_folder; // Unused variable;

/////////////////////////////
        // Run the simulation. 
        auto begin = std::chrono::high_resolution_clock::now(); // Timer.

        // Run model must take the json file as an input as well as the scenario name (scenario name is the quarantine input folder)
        std::vector<individual> residents = run_model(sim_number, scenario_name, parameters_json);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin); // Timer. 
    
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        
        // Filename to write too. 
        std::string filename = directory  + "/sim_number_" + std::to_string(sim_number) + ".csv";
    
        // Write output to file.
        std::ofstream output_file(filename);
        if(output_file.is_open()){

            output_file << "Individual, Household, Age, Current Vaccine, Current doses, Time of first dose, Time of last dose, Vaccine at infection, Doses at infection, Time latest dose at infection, Severity, Symptomatic, Time isolated, Detected case, Time of detection, Time of exposure, Time of infection, Time of symptom onset, Time of recovery, Cluster, Secondary Infections" << std::endl;
            
            int ind_num = 0;
            // To do: we must define both the vaccination status of the individual at infection and at end?
            for(individual & person: residents){
                bool is_infected = !std::isnan(person.covid.time_of_exposure); // Are they infected.
                vaccine_type  vac = person.vaccine_status.get_type();
                
                // Infected write!
                if(is_infected){
                    // Determine vaccine at time of infection.
                    vaccine_type infection_vac = person.infection_statistics.vaccine_status;
                    
                    output_file << ind_num << ", " << person.home_id << ", " << person.age << ", ";
                    
                    output_file << vac << ", " << person.vaccine_status.get_dose() <<", " << person.vaccine_status.get_first_time() << ", " << ((person.vaccine_status.get_dose()==2)?person.vaccine_status.get_time_of_vaccination():std::nan("1")) << ", "; // Vaccination at end of simulation.
                    
                    output_file << infection_vac << ", " << person.infection_statistics.doses << ", " << person.infection_statistics.time_of_last_dose << ", ";
                    
                    output_file << (person.covid.severe?"Severe":"Mild") << ", " << (person.covid.asymptomatic?"Asymptomatic":"Symptomatic") << ", "<< person.time_isolated << ", " << (person.covid.detected_case?"Detected":"Undetected") << ", " << person.covid.time_of_detection << ", ";
                    
                    output_file << person.covid.time_of_exposure << ", " << person.covid.time_of_infection << ", " << person.covid.time_of_symptom_onset << ", ";
                    output_file << person.covid.time_of_recovery << ", " << person.covid.cluster_number << ", " << person.who_infected.size() << std::endl;
                }
                ind_num ++;
            }
        // Close file.
        output_file.close();
        }
/////////////////////////////    
    
    return 0; 
}