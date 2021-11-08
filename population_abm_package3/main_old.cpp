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

void load_linelist(std::vector<double> & total_time, std::vector<double> & exposure_days, std::vector<double>& FOI, std::string filename){
    
    std::ifstream breach_linelist(filename);
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
}

class pseudo_individual{
    public:
    pseudo_individual(double t_start, double t_E, double FOI, int age_bracket_in, std::vector<double> & age_stratified_contacts,std::vector<double> & xi_k,int cluster_ref):start_exposure(t_start),finish_exposure(t_start + t_E),age_bracket(age_bracket_in),cluster(cluster_ref){
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
    double start_exposure;
    double finish_exposure;
    double alpha; // Scale to ensure expected secondary infections is FOI from quarantine. 
    int age_bracket; 
    int cluster; 
};

std::vector<individual> run_model(double beta_H, double beta_C, int num_houses, double average_house_size, std::vector<double> & age_brackets, std::vector<double> & population_pi, std::vector<std::vector<double>> &contact_matrix){
    
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
    double t_end = 300.0; // When do we finish? 
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

    

    // Quarantine breach load 
    std::vector<double> total_time;
    std::vector<double> FOI; 
    std::vector<double> exposure_days;
    load_linelist(total_time,exposure_days,FOI,"scenario1/arrivals_linelist_scenario1_sim1_test.csv");

    // Quarantine breach.
    for(int breach = 0; breach < total_time.size(); breach++){
        ++cluster_ref; 
        int age_bracket = rand()%16;
        double foi = FOI[breach];
        double t_E = exposure_days[breach];
        double time = total_time[breach];
        active_breaches.push_back(pseudo_individual(time, t_E, foi, age_bracket, contact_matrix[age_bracket], xi_k, cluster_ref));
        
    }
    
    for(auto breach: active_breaches){
    std::cout << breach.start_exposure << " ";
    std::cout << breach.finish_exposure << " ";
    std::cout << breach.cluster << " ";
    std::cout << breach.alpha << " active breaches \n";
    
    }

    std::vector<int> total_exposed_by_breach; total_exposed_by_breach.reserve(1000); // Magic number to minimise reallocation at first.
    std::cout << active_breaches[0].age_bracket << " " << active_breaches[0].alpha << " " << active_breaches[0].finish_exposure; 
    while(t < t_end){
        
        // Infection model!
        std::cout << "Current time = " << t << "\n";


        std::vector<int> exposed_by_breach; exposed_by_breach.reserve(1000); // Magic number to minimise reallocation at first.

        auto remove_breach = std::remove_if(active_breaches.begin(),active_breaches.end(),[&](auto breach)->bool{
            // Run a simulation of who the pseudo individual will contact. 
            // Who do they contact, will they become sick. 
            // The exposed individuals here... should not become exposed immediately, this time step they are not going to contribute to the FOI on the rest of the population.
            if(t >= breach.finish_exposure){ 
                return true; // Remove if their exposure period is over. 

            }else if(t >=breach.start_exposure && t < breach.finish_exposure){
                std::cout << t << std::endl;
            // return t >= breach.finish_exposure; // Remove them from active breaches if their exposure period is over. 
            auto & contacts = contact_matrix[breach.age_bracket];
            double alpha = breach.alpha;
            for(int i = 0; i < contacts.size(); i++){
                double Lambda = contacts[i];
                std::poisson_distribution<int> gen_num_contacts(alpha*Lambda*dt);

                size_t num_in_strata = age_matrix[i].size(); // Must be obtained from the reference to the age matrix.

                if(num_in_strata == 0){
                    continue;    // Skip strata as no-one is in it.
                }
        
                std::uniform_int_distribution<size_t> gen_reference(0,num_in_strata-1);
                int num_contacts = gen_num_contacts(generator);
                
                for(int j = 0; j < num_contacts; j++){
                    int contact_ref = gen_reference(generator);

                    individual & contact = residents[contact_ref]; // Define which community member was contacted.
                    // Track all contacts (We can then implement how good contact tracing is currently disabled)
                    if(contact.covid.infection_status == 'S'){
                    // Do they get infected.
                        double r = genunf_std(generator);

                        if(r < xi_k[i]){   // This part doesnt need the dt, because the dt is taken into account earlier by limiting the number of contacts per timesteps.
                            exposed_by_breach.push_back(contact_ref);
                            covid.seed_exposure(contact,t); // Random resident has become infected
                            contact.covid.cluster_number = breach.cluster;
                            
                        }
                    }   
                }
            }
        }
            return false;
 });

        active_breaches.erase(remove_breach, active_breaches.end());

        // Reinitialise the new symptomatic infections.
        std::vector<size_t> newly_symptomatic; newly_symptomatic.reserve(3000); // Reserve size so that reallocation isnt neccesary. Its a magic number.
        // t  +=dt;
        t = covid.covid_ascm(residents,houses,age_matrix,t,t+dt,dt,E_ref,I_ref,newly_symptomatic);
        
        // Now add the exposed by breach individuals. 
        E_ref.insert(E_ref.end(),exposed_by_breach.begin(),exposed_by_breach.end());

        std::cout << E_ref.size() << " exposed " << I_ref.size() << " infected \n"; 
        
    }

    return residents; // How many were exposed. 
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

    
    
    // Read in the scenario  files.
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
    double  scale_facility = parameters_json["scale_facility"];
    double  breach_rate = parameters_json["breach_rate"];

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
        std::vector<individual> residents = run_model(beta, beta_C, num_houses, house_size, age_brackets, population_pi, contact_matrix);
        
        // count_infections += output;
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
       
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        
        // Filename.
        // Increment by 1 to not start at 0, then adjust based on the "start_sim" input.
        int sim_number = i+1 + (start_sim-1);
        std::string filename = directory  + "/sim_number_" + std::to_string(sim_number) + ".csv";


    // std::cout << static_cast<double>(count_infections)/ static_cast<double>(num_sims) << " average secondary infections \n";
    // Write file.
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

        };

    }


    return 0;
 
};



