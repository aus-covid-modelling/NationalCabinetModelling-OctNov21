#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include "../vaccine_abm/individual.h"
#include "../vaccine_abm/households.h"
#include "../vaccine_abm/ibm_simulation.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include <algorithm>
#include "../nlohmann/json.hpp"
#include "../vaccine_abm/USER_random.h"



class pseudo_individual{
    public:
    pseudo_individual(double t, double t_E, double FOI, int age_bracket_in, std::vector<double> & age_stratified_contacts,std::vector<double> & xi_k):finish_exposure(t + t_E),age_bracket(age_bracket_in){
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
    double finish_exposure;
    double alpha; // Scale to ensure expected secondary infections is FOI from quarantine. 
    int age_bracket; 
};

std::pair<std::vector<double>,std::vector<double>> run_model(double beta_H, double beta_C, int num_houses, double average_house_size, std::vector<double> & age_brackets, std::vector<double> & population_pi, std::vector<std::vector<double>> &contact_matrix){
    
    // Number of age brackets.
    size_t num_brackets = population_pi.size();
    
    // Disease model initialisation 
    std::vector<double> b{100.0,101.0};// Testing with no ttiq + 100 days.  
    std::vector<double> w{1.0}; // Just to make sure it works. 
    disease_model covid(beta_H, beta_C, contact_matrix,b,w); // Load the disease model
    
    // Storage for households and residents.
    std::vector<household> houses; houses.reserve(num_houses);
    std::vector<individual> residents; residents.reserve(std::floor(num_houses*average_house_size));

    // Active quarantine breaches. 
    std::vector<pseudo_individual> active_breaches; 

    // Read in the Vaccine parameters files.
    std::ifstream vaccine_in("vaccine_parameters.json");
    nlohmann::json vaccine_json;
    vaccine_in >> vaccine_json;
    vaccine_in.close(); // Close the file.

    std::vector<double> xi_k = vaccine_json["Unvaccinated"]["xi"]; // This is for the pseudo individual. 

    vaccine_parameters no_vaccine(0,vaccine_type::none, vaccine_json["Unvaccinated"]["xi"], vaccine_json["Unvaccinated"]["tau"], vaccine_json["Unvaccinated"]["symptomatic"]); // No vaccination will be passed into everyone. Required for the start of the simulation.
    

    // Initialise all of the community...
    community city(0,num_houses,average_house_size,residents,houses,age_brackets,population_pi,no_vaccine); // Vaccinated proportion not used here now.
    std::cout << "residents = " << residents.size() << " houses = " << houses.size() << "\n";
   
    // Timestepping parameters
    double t = 0.0;
    double t_end = 30.0; // When do we finish? 
    double dt = pow(2.0,-2.0); // Probably doesnt need 16 per day but it runs so quick...
    
    // Create memory that tracks who is exposed, E_ref, and who is infected, I_ref. gen_res is used to sample from the list of residents uniformly.
    std::vector<size_t> E_ref; E_ref.reserve(10000); // Magic number reserving memory.
    std::vector<size_t> I_ref; I_ref.reserve(10000); // Magic number of reserved.
    
    // Assemble the age_matrix (this is a list of people that are in each age bracket).
    std::vector<std::vector<int>> age_matrix(num_brackets);
    assemble_age_matrix(residents,age_matrix); // Nobody moves from the age matrix so only have to do it once.
    
    std::uniform_int_distribution<size_t> gen_res(0,residents.size()-1);

    int cluster_ref = 0; // Track the exposure number, can show that one dominates. Do not update cluster ref here. This is a reference for original infections. 
    // int initial_infections = 0; // Count initial infections. 
    // while(initial_infections < 1){

    //     int exposed_resident = gen_res(generator); // Randomly sample from all the population.

    //         if(residents[exposed_resident].vaccine_status.get_type()==vaccine_type::none){
    //             if(residents[exposed_resident].covid.infection_status!='E'){
    //                 covid.seed_exposure(residents[exposed_resident],t); // Random resident has become infected
    //                 residents[exposed_resident].covid.cluster_number = cluster_ref;
    //                 ++initial_infections;
    //                 E_ref.push_back(exposed_resident); // Start tracking them.
    //             }
    //         }
    // }

    // For checking secondary infections Tost and GI. 
    int exposed_resident = gen_res(generator); // Randomly sample from all the population.
    if(residents[exposed_resident].covid.infection_status!='E'){
        covid.seed_exposure(residents[exposed_resident],t); // Random resident has become infected
        residents[exposed_resident].covid.cluster_number = cluster_ref;
        E_ref.push_back(exposed_resident); // Start tracking them.
    }
    

    // Quarantine breach check. 
    // Have just one quarantine breach for testing. 
    int age_bracket = 1; // Sample age bracket. 
    double FOI = 1.0;
    double t_E = 10.0;
    active_breaches.push_back(pseudo_individual(t, t_E, FOI, age_bracket, contact_matrix[age_bracket], xi_k)); 

    while(t < t_end){

        // Infection model!
        std::cout << "Current time = " << t << "\n";

            // Quarantine breach check. 
        bool traveller_breach = false; // Disable for testing
        if(traveller_breach){
            // Get a value from the line list.
            int age_bracket = 0; // Sample age bracket. 
            double FOI = 1.0;
            double t_E = 10.0;
            active_breaches.push_back(pseudo_individual(t, t_E, FOI, age_bracket, contact_matrix[age_bracket], xi_k)); 
            
        }
        bool worker_breach = false; // Disable for testing
        if(worker_breach){
            // Get a value from the line list.
            int age_bracket = 1; // Sample age bracket. 
            double FOI = 1.0;
            double t_E = 10.0;
            active_breaches.push_back(pseudo_individual(t, t_E, FOI, age_bracket, contact_matrix[age_bracket], xi_k)); 
            
        }

        std::vector<int> exposed_by_breach; exposed_by_breach.reserve(1000); // Magic number to minimise reallocation at first. 
        std::remove_if(active_breaches.begin(),active_breaches.end(),[&](auto breach)->bool{
            // Run a simulation of who the pseudo individual will contact. 
            // Who do they contact, will they become sick. 
            // The exposed individuals here... should not become exposed immediately, this time step they are not going to contribute to the FOI on the rest of the population. 

            return t >= breach.finish_exposure; // Remove them from active breaches if their exposure period is over. 
        });

        // Reinitialise the new symptomatic infections.
        std::vector<size_t> newly_symptomatic; newly_symptomatic.reserve(3000); // Reserve size so that reallocation isnt neccesary. Its a magic number.
        
        t = covid.covid_ascm_R0(residents,houses,age_matrix,t,t+dt,dt,E_ref,I_ref,newly_symptomatic);
        
        // Now add the exposed by breach individuals. 
        E_ref.insert(E_ref.end(),exposed_by_breach.begin(),exposed_by_breach.end());

        std::cout << E_ref.size() << " exposed " << I_ref.size() << " infected \n"; 
        
    }

    // Finished model.
    std::vector<double>tost;
    std::vector<double> gi;
    double time_of_symptoms = residents[exposed_resident].covid.time_of_symptom_onset;
    // double time_of_infection = residents[exposed_resident].covid.time_of_infection;
    double time_of_infection = residents[exposed_resident].covid.time_of_exposure;
    for(auto & infectee: residents[exposed_resident].who_infected){
        // gi.push_back(residents[infectee].covid.time_of_infection - time_of_infection);
        gi.push_back(residents[infectee].covid.time_of_exposure - time_of_infection);
        tost.push_back(residents[infectee].covid.time_of_exposure - time_of_symptoms);
    }
    // return residents[exposed_resident].who_infected.size();
    return make_pair(tost,gi);
}

// The main script driver that runs and call run_model.
int main(int argc, char *argv[]){

    if(argc == 1){
        std::cout << "ERROR: Parameter file not loaded!" << std::endl;
        return 0;
    }
    // Parameter read in and simulation here!
    std::string filename(argv[1]); // Reads in json parameter. 
    int start_sim;

    if(argc < 3){throw std::logic_error("Sim number not provided.\n");}
    
    // Read sim number
    std::istringstream iss( argv[2] );
    if(iss >> start_sim)
    {
        std::cout << "Sim number " << start_sim << std::endl;
    }else{
        throw std::logic_error("Sim number failed to convert\n");
    }

    
    
    // Read in the parameter files.
    std::ifstream parameters_in(filename);


    nlohmann::json parameters_json;
    parameters_in >> parameters_json;
    parameters_in.close(); // Close the file.

    // Assign values from json! 
    int num_sims = parameters_json["num_sims"];
    std::cout << "Number of simulations: " << num_sims << std::endl;

    // Load paramters. 
    double  beta = parameters_json["Beta_H"];
    double  TP = parameters_json["TP"];
    double  house_size = parameters_json["average_house_size"];
    int     num_houses = parameters_json["num_houses"];

    std::string ttiq_type = parameters_json["ttiq"];
    std::string vaccine_filenames = parameters_json["vaccine_inputs"];
    
    // Contact matrix filename. 
    std::string contact_matrix_filename = parameters_json["contact_matrix"];

    // // Create folder.
    std::string folder = parameters_json["folder_name"];
    std::string directory = (std::string) parameters_json["output_directory"] + folder;

    #ifdef _WIN32
        int top_folder = mkdir(((std::string) parameters_json["output_directory"]).c_str());
        int main_folder = mkdir(directory.c_str()); // Create folder.
    #else
        int top_folder = mkdir(((std::string) parameters_json["output_directory"]).c_str(), 0777);
        int main_folder = mkdir(directory.c_str(),0777); // Create folder.
    #endif
    (void) main_folder; // Unused variable;
    

    //  Population demographic parameters (There should be 16 here)
    // Now found in the parameters file (field population_pi)
    std::vector<double> population_pi = parameters_json["population_pi"];
    int num_brackets    = (int) population_pi.size();
    
    nlohmann::json & vaccinated_proportion = parameters_json["vaccinated_proportion"];
    std::vector<double> pfizer_proportion = vaccinated_proportion["Pfizer"];
    std::vector<double> az_proportion = vaccinated_proportion["Astrazeneca"];
    std::vector<double> moderna_proportion = vaccinated_proportion["Moderna"];
    
    for(auto x : pfizer_proportion) std::cout << x << ", ";
    std::cout << "\n ";for(auto x : az_proportion) std::cout << x << ", ";
    std::cout << "\n ";for(auto x : moderna_proportion) std::cout << x << ", ";

    // Define the age brackets that are used in the contact_matrix (this is used to generate an age for each individual)
    std::vector<double> age_brackets{0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
    if(age_brackets.size()!=num_brackets){
        throw std::logic_error("Error: age brackets not the same size as population proportion!!");
    }
    
    // Contact matrix (currently importing Prem's Australia matrix).
    std::vector<std::vector<double>> contact_matrix(num_brackets,std::vector<double>(num_brackets,0.0));
    std::ifstream matrix_read;
    matrix_read.open(contact_matrix_filename);
    if(matrix_read.is_open()){
    for(size_t i = 0; i < num_brackets;i++){
        for(size_t j = 0; j <num_brackets;j++){
            matrix_read >> contact_matrix[i][j];
        }
    }
    matrix_read.close();
    }else{
        throw std::logic_error(contact_matrix_filename + " not found in working directory.\n");
    }
    
    
    // The R0 calculation that has been updated to account for all of the heterogeneity built in at the moment....
        // Read in the Vaccine parameters files.
    std::ifstream vaccine_in("vaccine_parameters.json");
    nlohmann::json vaccine_json;
    vaccine_in >> vaccine_json;
    vaccine_in.close(); // Close the file.

    double tau_s = vaccine_json["Unvaccinated"]["tau"][0]; double tau_sc = vaccine_json["Unvaccinated"]["tau"][1];

    vaccine_parameters no_vaccine(0,vaccine_type::none, vaccine_json["Unvaccinated"]["xi"], vaccine_json["Unvaccinated"]["tau"], vaccine_json["Unvaccinated"]["symptomatic"]); // No vaccination will be passed into everyone. Required for the start of the simulation.

    double sum_expression = 0.0;
    for(int k = 0; k < (int) num_brackets;k++){
        double xi_k = no_vaccine.get_susceptibility(k);
        
        double internal_sum = 0.0;
        for(int i = 0; i < (int) num_brackets; i++){
            double lambda_ik = contact_matrix[i][k];
            internal_sum += lambda_ik*((tau_s - tau_sc)*no_vaccine.get_proportion_asymptomatic(i) + tau_sc)*population_pi[i];
        };
        sum_expression += internal_sum*xi_k;
    };
    

    std::cout << num_sims << " simulations " << std::endl;
    int count_infections = 0;
    std::vector<double> TOTAL_TOST;
    std::vector<int> total_si;
    std::vector<double> TOTAL_GI;
    for(int i = 0; i < num_sims; i++){
        // Calculate beta_C for the simulation.
        std::cout << "Sim Number = " << i+1 << "\n";
        auto begin = std::chrono::high_resolution_clock::now();
        // double beta_C = TP/((5.0)*sum_expression);
        double beta_C = TP/(sum_expression*(2.29*(5.1-2.5) + 1.5));

        std::cout << beta_C << " beta_C IS HARDCODED BE CAREFUL \n";
        std::pair<std::vector<double>,std::vector<double>> output = run_model(beta, beta_C, num_houses, house_size, age_brackets, population_pi, contact_matrix);
        std::vector<double> tost = output.first;
        std::vector<double> gi = output.second;
        count_infections+=tost.size();
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
        TOTAL_TOST.insert(TOTAL_TOST.end(),tost.begin(),tost.end()); 
        TOTAL_GI.insert(TOTAL_GI.end(),gi.begin(),gi.end()); 
        total_si.push_back(tost.size());
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        
        // Filename.
        // Increment by 1 to not start at 0, then adjust based on the "start_sim" input.
        int sim_number = i+1 + (start_sim-1);
        std::string filename = directory  + "/sim_number_" + std::to_string(sim_number) + ".csv";


        // // Write file.
        // std::ofstream output_file(filename);
        // if(output_file.is_open()){

        //     output_file << "Individual, Household, Age, Current Vaccine, Current doses, Time of first dose, Time of last dose, Vaccine at infection, Doses at infection, Time latest dose at infection, Severity, Symptomatic, Time isolated, Detected case, Time of detection, Time of exposure, Time of infection, Time of symptom onset, Time of recovery, Cluster, Secondary Infections" << std::endl;
            
        //     int ind_num = 0;
        //     // To do: we must define both the vaccination status of the individual at infection and at end?
        //     for(individual & person: residents){
        //         bool is_infected = !std::isnan(person.covid.time_of_exposure); // Are they infected.
        //         vaccine_type  vac = person.vaccine_status.get_type();
                
        //         // Infected write!
        //         if(is_infected){
        //             // Determine vaccine at time of infection.
        //             vaccine_type infection_vac = person.infection_statistics.vaccine_status;
                    
        //             output_file << ind_num << ", " << person.home_id << ", " << person.age << ", ";
                    
        //             output_file << vac << ", " << person.vaccine_status.get_dose() <<", " << person.vaccine_status.get_first_time() << ", " << ((person.vaccine_status.get_dose()==2)?person.vaccine_status.get_time_of_vaccination():std::nan("1")) << ", "; // Vaccination at end of simulation.
                    
        //             output_file << infection_vac << ", " << person.infection_statistics.doses << ", " << person.infection_statistics.time_of_last_dose << ", ";
                    
        //             output_file << (person.covid.severe?"Severe":"Mild") << ", " << (person.covid.asymptomatic?"Asymptomatic":"Symptomatic") << ", "<< person.time_isolated << ", " << (person.covid.detected_case?"Detected":"Undetected") << ", " << person.covid.time_of_detection << ", ";
                    
        //             output_file << person.covid.time_of_exposure << ", " << person.covid.time_of_infection << ", " << person.covid.time_of_symptom_onset << ", ";
        //             output_file << person.covid.time_of_recovery << ", " << person.covid.cluster_number << ", " << person.who_infected.size() << std::endl;
        //         }
        //         ind_num ++;
        //     }


        // // Close file.
        //  output_file.close();

        // };

    }

          std::ofstream output_file("TOST_1_5.csv");
        if(output_file.is_open()){
            output_file << "TOST\n";
            for(auto x: TOTAL_TOST){
                output_file << x << "\n";
            }
        }
        output_file.close();

        std::ofstream output_gi("GI_1_5.csv");
        if(output_gi.is_open()){
            output_gi << "GI\n";
            for(auto x: TOTAL_GI){
                output_gi << x << "\n";
            }
        }
        output_gi.close();

                  std::ofstream si("SI_1_5.csv");
        if(si.is_open()){
            si << "SI\n";
            for(auto x: total_si){
                si << x << "\n";
            }
        }
        si.close();

    std::cout << static_cast<double>(count_infections)/ static_cast<double>(num_sims) << " average secondary infections \n";
    
    return 0;
 
};



