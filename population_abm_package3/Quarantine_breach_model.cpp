#include "Quarantine_breach_model.h"
#include "Quarantine_outputs.h"
#include <fstream>
#include <sstream>
#include "../vaccine_abm/nbinrnd.h"

std::vector<individual> run_model(int sim_number,std::string scenario_name, nlohmann::json & inputs){
    
    // File locations for TTIQ and Vaccine parameters. 
    std::string ttiq_type = inputs["ttiq"];
    std::string vaccine_filenames = inputs["vaccine_inputs"];

    //  Population demographic parameters 
    std::vector<double> population_pi = inputs["population_pi"];    
    size_t num_brackets = population_pi.size();
    std::vector<double> age_brackets = inputs["age_brackets"];
    if(age_brackets.size()!=num_brackets){
        throw std::logic_error("Error: age brackets not the same size as population proportion!!");
    }

    // Contact matrix filename. 
    std::string contact_matrix_filename = inputs["contact_matrix"];

    // Contact matrix - read in. 
    std::vector<std::vector<double>> contact_matrix(num_brackets,std::vector<double>(num_brackets,0.0));
    std::ifstream matrix_read(contact_matrix_filename);
    if(matrix_read.is_open()){
        std::string line;
        double value;
        for(int i = 0; i < num_brackets; i++){
            std::getline(matrix_read,line); // Read the line. 
            std::stringstream stream_line(line);
            std::string row_val;
            for(int j = 0; j < num_brackets; j++){
                std::getline(stream_line,row_val,','); // Row val is the corresponding double we are after. 
                std::stringstream stream_row(row_val);
                stream_row >> value;
                contact_matrix[i][j] = value;
                // std::cout << value << ", ";
                
            }
            // std::cout << std::endl;
        }
        matrix_read.close();
    } else {
        throw std::logic_error("The contact matrix file " + contact_matrix_filename + " was not found.");
    }

     std::cout << "Contact matrix " << std::endl;
    for(int i = 0; i < contact_matrix.size(); i++){
        for(int j = 0;j < contact_matrix[i].size(); j++){
            std::cout << contact_matrix[i][j] << ", ";
        }
        std::cout << std::endl;
    }

    // Load the TTIQ distributions. 
    std::ifstream ttiq_read("ttiq_distributions.csv");
    std::vector<double> ttiq_days, partial, optimal;
    if(ttiq_read.is_open()){
        std::string line;
        double value;
        std::getline(ttiq_read,line); // Get first line. 
        // Do nothing this is a title.
        // Second line onwards. 
        while(std::getline(ttiq_read,line)){
            std::stringstream stream_line(line);
            std::string row_val;
            // std::vector<double> row;
            for(int col_number = 0; col_number < 3; col_number++){
                std::getline(stream_line,row_val,','); // Row val is the corresponding double we are after. 
                std::stringstream stream_row(row_val);
                stream_row >> value;
                if(col_number ==0){
                    ttiq_days.push_back(value);
                }
                if(col_number == 1){
                    partial.push_back(value);
                }
                if(col_number ==2){
                    optimal.push_back(value);
                }
            }
        }
        ttiq_read.close();
    } else {
        throw std::logic_error("The ttiq distributions file was not found.");
    }
    
    // Error check. 
    if(ttiq_days.size()!=partial.size() || ttiq_days.size()!= optimal.size()){
        throw std::logic_error("Dimension mismatch in ttiq distributions! \n");
    }

    

    // Read in the Vaccine parameters files.
    std::ifstream vaccine_in(vaccine_filenames);
    nlohmann::json vaccine_json;
    vaccine_in >> vaccine_json; // Assign the vaccine information into the model. 
    vaccine_in.close(); // Close the file.

    // Vaccine parameters being assigned. 
    vaccine_parameters no_vaccine(0,vaccine_type::none, vaccine_json["Unvaccinated"]["xi"], vaccine_json["Unvaccinated"]["tau"], vaccine_json["Unvaccinated"]["symptomatic"]); // No vaccination will be passed into everyone. Required for the start of the simulation.
    std::vector<double> xi_k = vaccine_json["Unvaccinated"]["xi"]; // This is for the pseudo individual constructor. 
    
    // Calculate the beta_C for a completely susceptible population from the TP. 
    double tau_s = vaccine_json["Unvaccinated"]["tau"][0]; double tau_sc = vaccine_json["Unvaccinated"]["tau"][1]; // Assign tau for book keeping.
    std::vector<double> alpha = inputs["relative_infectiousness"];
    

    double sum_expression = 0.0;
    for(int k = 0; k < (int) num_brackets;k++){
        double xi_k = no_vaccine.get_susceptibility(k);
        
        double internal_sum = 0.0;
        for(int i = 0; i < (int) num_brackets; i++){
            double lambda_ik = contact_matrix[i][k];
            internal_sum += alpha[i]*lambda_ik*((tau_s - tau_sc)*no_vaccine.get_proportion_asymptomatic(i) + tau_sc)*population_pi[i];
        };
        sum_expression += internal_sum*xi_k;
    };
    // Read other in inputs from json!
    double  beta_H = inputs["Beta_H"];
    double  TP = inputs["TP"];
    double  house_size = inputs["average_house_size"];
    int     num_houses = inputs["num_houses"];
    // Assign beta_C (used for calculating community infection rates).
    double beta_scale = TP/(sum_expression*(2.29*(5.1-2.5) + 1.5)); // This is hardcoded, be careful if anything changes. 
    std::vector<double>beta_C = alpha;
    for(auto & x:beta_C){ // Scale so that it is the appropriate size. 
        x = beta_scale*x;
    }

    // Disease model initialisation 
    std::vector<double> b(ttiq_days.begin(),ttiq_days.end());
    std::vector<double> w;
    if(ttiq_type=="partial"){
        w = std::vector<double>(partial.begin(),partial.end());
    }else if(ttiq_type == "optimal"){
        w = std::vector<double>(optimal.begin(), optimal.end());
    }else{
        throw std::logic_error("Unrecognised TTIQ_type, please choose either partial or optimal (CASE SENSITIVE)\n");
    }
    // Disease model.
    disease_model covid(std::vector<double>(num_brackets,beta_H), beta_C, contact_matrix, b, w); // Load the disease model

    // Storage for households and residents.
    std::vector<household> houses; houses.reserve(num_houses);
    std::vector<individual> residents; residents.reserve(std::floor(num_houses*house_size));

    // Active quarantine breaches. 
    std::vector<pseudo_individual> active_breaches; 

    // Initialise all of the community...
    community city(0,num_houses,house_size,residents,houses,age_brackets,population_pi,no_vaccine); // Vaccinated proportion not used here now.
    std::cout << "residents = " << residents.size() << " houses = " << houses.size() << "\n";
   
    // Vaccinated proportion. 
    nlohmann::json & vaccinated_proportion = inputs["vaccinated_proportion"];
    std::vector<double> pfizer_first = vaccinated_proportion["Pfizer"]["First_Dose"];
    std::vector<double> az_first = vaccinated_proportion["Astrazeneca"]["First_Dose"];
    std::vector<double> moderna_first = vaccinated_proportion["Moderna"]["First_Dose"];
    std::vector<double> pfizer_second = vaccinated_proportion["Pfizer"]["Second_Dose"];
    std::vector<double> az_second = vaccinated_proportion["Astrazeneca"]["Second_Dose"];
    std::vector<double> moderna_second = vaccinated_proportion["Moderna"]["Second_Dose"];

    std::cout << std::endl;
    
    nlohmann::json & pfizer_input = vaccine_json["Pfizer"];
    nlohmann::json & moderna_input = vaccine_json["Moderna"];
    nlohmann::json & astrazeneca_input = vaccine_json["Astrazeneca"];

    // Construct vaccination parameters 
    std::vector<vaccine_parameters> pfizer{vaccine_parameters(1,vaccine_type::pfizer,pfizer_input["First_Dose"]["xi"],pfizer_input["First_Dose"]["tau"],pfizer_input["First_Dose"]["symptomatic"]),vaccine_parameters(2,vaccine_type::pfizer,pfizer_input["Second_Dose"]["xi"],pfizer_input["Second_Dose"]["tau"],pfizer_input["Second_Dose"]["symptomatic"])};

    std::vector<vaccine_parameters> moderna{vaccine_parameters(1,vaccine_type::moderna,moderna_input["First_Dose"]["xi"],moderna_input["First_Dose"]["tau"],moderna_input["First_Dose"]["symptomatic"]),vaccine_parameters(2,vaccine_type::moderna,moderna_input["Second_Dose"]["xi"],moderna_input["Second_Dose"]["tau"],moderna_input["Second_Dose"]["symptomatic"])};

    std::vector<vaccine_parameters> astrazeneca{vaccine_parameters(1,vaccine_type::astrazeneca,astrazeneca_input["First_Dose"]["xi"],astrazeneca_input["First_Dose"]["tau"],astrazeneca_input["First_Dose"]["symptomatic"]),vaccine_parameters(2,vaccine_type::astrazeneca,astrazeneca_input["Second_Dose"]["xi"],astrazeneca_input["Second_Dose"]["tau"],astrazeneca_input["Second_Dose"]["symptomatic"])};
     
    // Vaccinate the individuals - Will the proportion vaccinated have to be different
    for(int i = 0; i < residents.size(); i++){// For each individual.
        size_t dim_age_band; // Quantium related. 
        double age = residents[i].age;
        double pre_time = -20;
        if(age < 12){
            continue; // Skip the individual
        }else if(age >= 12 && age < 16){
            dim_age_band = 0;
        }else if(age >=80){
            dim_age_band = 8;
        }else{
            dim_age_band = std::floor(residents[i].age/10.0);    // Age bracket.
        }

        double r = genunf_std(generator);
        // Check which vaccine they want - based off the inputs.  
        if(r < pfizer_first[dim_age_band]){
            // Pfizer first dose.
            residents[i].vaccine_status.vaccinate_individual(pfizer[0], dim_age_band, pre_time);

        }else if(r < pfizer_first[dim_age_band] + pfizer_second[dim_age_band]){
            // Pfizer second dose.
            residents[i].vaccine_status.vaccinate_individual(pfizer[0], dim_age_band, pre_time);
            residents[i].vaccine_status.vaccinate_individual(pfizer[1], dim_age_band, pre_time);

        }else if(r < pfizer_first[dim_age_band] + pfizer_second[dim_age_band] + moderna_first[dim_age_band]){
            // Moderna First dose.
            residents[i].vaccine_status.vaccinate_individual(moderna[0], dim_age_band, pre_time);

        }else if(r < pfizer_first[dim_age_band] + pfizer_second[dim_age_band] + moderna_first[dim_age_band] + moderna_second[dim_age_band]){
            // Moderna second dose. 
            residents[i].vaccine_status.vaccinate_individual(moderna[0], dim_age_band, pre_time);
            residents[i].vaccine_status.vaccinate_individual(moderna[1], dim_age_band, pre_time);

        }else if(r < pfizer_first[dim_age_band] + pfizer_second[dim_age_band] + moderna_first[dim_age_band] + moderna_second[dim_age_band] + az_first[dim_age_band]){
            // az First dose. 
            residents[i].vaccine_status.vaccinate_individual(astrazeneca[0], dim_age_band, pre_time);

        }else if(r < pfizer_first[dim_age_band] + pfizer_second[dim_age_band] + moderna_first[dim_age_band] + moderna_second[dim_age_band] + az_first[dim_age_band] + az_second[dim_age_band]){
            // AZ second dose. 
            residents[i].vaccine_status.vaccinate_individual(astrazeneca[0], dim_age_band, pre_time);
            residents[i].vaccine_status.vaccinate_individual(astrazeneca[1], dim_age_band, pre_time);

        }else{
            // Do nothing - you are unvaccinated!
        }
    }

    // Error check!
    for(int dim_age_band=0; dim_age_band<9; dim_age_band++){
        if(pfizer_first[dim_age_band] + pfizer_second[dim_age_band] + moderna_first[dim_age_band] + moderna_second[dim_age_band] + az_first[dim_age_band] + az_second[dim_age_band]>1){
            throw std::logic_error("The probability of being vaccinated is greater than 1? ... Check proportions vaccinated! \n");
        }
    }
    
    // Timestepping parameters
    double t = 0.0;
    double t_end = inputs["end_time"]; // When do we finish? 
    double dt = pow(2.0,-2.0); // Probably doesnt need 16 per day but it runs so quick...
    
    // Create memory that tracks who is exposed, E_ref, and who is infected, I_ref. gen_res is used to sample from the list of residents uniformly.
    std::vector<size_t> E_ref; E_ref.reserve(10000); // Magic number reserving memory.
    std::vector<size_t> I_ref; I_ref.reserve(10000); // Magic number of reserved.
    
    // Assemble the age_matrix (this is a list of people that are in each age bracket).
    std::vector<std::vector<int>> age_matrix(num_brackets);
    assemble_age_matrix(residents,age_matrix); // Nobody moves from the age matrix so only have to do it once.

    // Use cluster ref to track the infections phylogenetic tree.
    std::uniform_int_distribution<size_t> gen_res(0,residents.size()-1); 
    int cluster_ref = 0; 
    int initial_infections = 0; // Count initial infections.
    int total_initial_infected = inputs["initial_infections"]; 
    while(initial_infections < total_initial_infected){
        int exposed_resident = gen_res(generator); // Randomly sample from all the population.
        if(residents[exposed_resident].vaccine_status.get_type()==vaccine_type::none){
            if(residents[exposed_resident].covid.infection_status!='E'){
                covid.seed_exposure(residents[exposed_resident],t); // Random resident has become infected
                residents[exposed_resident].covid.cluster_number = cluster_ref;
                ++initial_infections;
                E_ref.push_back(exposed_resident); // Start tracking them.
            }
        }
    }

    // Quarantine breach age distribution. 
    std::vector<double> adult_age_distribution = inputs["adult_age_distribution"];
    std::vector<double> kid_age_distribution = inputs["kid_age_distribution"];
    std::discrete_distribution<int>gen_adult_age(adult_age_distribution.begin(),adult_age_distribution.end()); 
    std::discrete_distribution<int>gen_child_age(kid_age_distribution.begin(),kid_age_distribution.end()); 

    // Load Quarantine breach linelist 
    std::vector<double> total_time;
    std::vector<double> FOI; 
    std::vector<double> exposure_days;
    std::vector<int> kid;
    load_linelist(total_time,exposure_days,FOI,kid, "arrivals_linelist/" + scenario_name + "/arrivals_linelist_" + scenario_name + "_sim" + std::to_string(sim_number) + ".csv");

    // Store the quarantine breaches in a vector of active breaches. 
    double time_shift = inputs["time_shift"]; // This shifts the time of the breach so that we have the appropriate case numbers at the starting point. 
    for(int breach = 0; breach < total_time.size(); breach++){
        ++cluster_ref; 
        double foi = FOI[breach];
        double t_E = exposure_days[breach];
        double time = total_time[breach] + time_shift;
        int age_bracket;
    
        if(kid[breach]==1){
            age_bracket = gen_child_age(generator);
            // foi = children_lower_infectious*foi;
        }else{
            age_bracket = gen_adult_age(generator);
        }

        // if(kid[breach]==1){
        // std::cout << t_E << " " << time << " " << age_bracket << " "<< kid[breach] << " " << foi << std::endl;
        // }
        

        active_breaches.push_back(pseudo_individual(time, t_E, foi, age_bracket, contact_matrix[age_bracket], xi_k, cluster_ref));   
    }


    while(t < t_end){
        
        // Infection model!
        std::cout << "Current time = " << t << "\n";

        // Simulate all the breach events. 
        std::vector<int> exposed_by_breach; exposed_by_breach.reserve(100); // Magic number to minimise reallocation at first.
        auto remove_breach = std::remove_if(active_breaches.begin(),active_breaches.end(),[&](auto breach)->bool{
            // Run a simulation of who the pseudo individual will contact. 
            // Who do they contact, will they become sick. 
            // The exposed individuals here... should not become exposed immediately, this time step they are not going to contribute to the FOI on the rest of the population.
            if(t >= breach.finish_exposure){ 
                return true; // Remove if their exposure period is over. 
            }else if(t >=breach.start_exposure && t < breach.finish_exposure){
                // return t >= breach.finish_exposure; // Remove them from active breaches if their exposure period is over. 
                auto & contacts = contact_matrix[breach.age_bracket];
                double alpha = breach.alpha;

                for(int i = 0; i < contacts.size(); i++){
                    double Lambda = contacts[i];

                    // Update for negative binomial. 
                    double r_daily = 0.050;
                    double mu_contacts = alpha*Lambda; //
                    double p = mu_contacts/(r_daily+mu_contacts);

                    // Note that contact matrix is square so doesnt matter which dimension has their size checked. 
                    nbinrnd gen_num_contacts(r_daily*dt,p); // Define the Distribution from which we sample community contacts.



                    // std::poisson_distribution<int> gen_num_contacts(alpha*Lambda*dt);
                    size_t num_in_strata = age_matrix[i].size(); // Must be obtained from the reference to the age matrix.

                    if(num_in_strata == 0){
                        continue;    // Skip strata as no-one is in it.
                    }
            
                    std::uniform_int_distribution<size_t> gen_reference(0,num_in_strata-1);
                    int num_contacts = gen_num_contacts(generator);
                    
                    for(int j = 0; j < num_contacts; j++){
                        int contact_ref = gen_reference(generator);

                        individual & contact = residents[contact_ref]; // Define which community member was contacted.
                        
                        if(contact.covid.infection_status == 'S'){
                        // Do they get infected.
                            double r = genunf_std(generator);

                            if(r < xi_k[i]){   // This part doesnt need the or the alpha dt, because the dt is taken into account earlier by limiting the number of contacts per timesteps. xi is always less than 1. 
                                exposed_by_breach.push_back(contact_ref); // List of individuals exposed by breach. 
                                covid.seed_exposure(contact,t); // Random resident has become infected
                                contact.covid.cluster_number = breach.cluster; // Pass on the cluster number so we can determine exactly who was infected by each breach event. 
                            }
                        }   
                    }
                }
            }
            return false;
        });
        active_breaches.erase(remove_breach, active_breaches.end()); // Delete any breaches that have finished. 

        // Reinitialise the newly symptomatic infections - set to empty.
        std::vector<size_t> newly_symptomatic; newly_symptomatic.reserve(3000); // Reserve size so that reallocation isnt always neccesary. Its a magic number.

        t = covid.covid_ascm(residents,houses,age_matrix,t,t+dt,dt,E_ref,I_ref,newly_symptomatic); // Do not change the t + dt here, otherwise it will break the relationship between breach events and covid simulation. 
        
        // Now add the exposed by breach individuals to the exposed reference vector. 
        E_ref.insert(E_ref.end(),exposed_by_breach.begin(),exposed_by_breach.end());
        std::cout << E_ref.size() << " exposed " << I_ref.size() << " infected \n"; 
        
    }

    return residents; // Return residents (move constructor) to write to file. 
}
