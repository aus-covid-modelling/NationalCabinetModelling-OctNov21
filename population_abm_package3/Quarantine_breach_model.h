#ifndef QUARANTINE_BREACH_MODEL
#define QUARANTINE_BREACH_MODEL
#include <vector>
#include "../nlohmann/json.hpp"
#include "../vaccine_abm/individual.h"
#include "../vaccine_abm/households.h"
#include "../vaccine_abm/ibm_simulation.h"
#include "../vaccine_abm/USER_random.h"
#include "../vaccine_abm/vaccine.h"
std::vector<individual> run_model(int sim_number, std::string scenario_name, nlohmann::json & inputs);
#endif