function AgentCharacteristics = initialise_infection_status(AgentCharacteristics,parameters)
 
% This script initialises infection and vaccination status of all agents

% AgentCharacteristics.InfectionStatus(i,1)= 
% 0: susceptible,
% 1: exposed, 
% 2: presymptomatic infectious, 
% 3: symptomatic infectious, 
% 4: presymptomatic asym infectious, 
% 5: asym infectious, 
% 6: recovered. 

% AgentCharacteristics.InfectionStatus(i,2)=
% Time to leave an E compartment 
% AgentCharacteristics.InfectionStatus(i,3)=
% Time to leave an PI/PA compartment
% AgentCharacteristics.InfectionStatus(i,4)=
% Time to leave an I/A compartment

AgentCharacteristics.InfectionStatus = zeros(parameters.PopSize, 4);

% Current test sensitivity of agents
% AgentCharacteristics.TestSensitivity(i,1)= 
% 1: if probability of true positive test for agent i is >0, 
% 0: otherwise
AgentCharacteristics.TestSensitivity = zeros(parameters.PopSize, 1);

% Contacts{j,1} = [ (unique agent IDs of contacts of agent j while infectious) binary_value timestep]
% Where the binary_value = 1 if it is a household contact, 0 otherwise
AgentCharacteristics.Contacts = cell(parameters.PopSize,1);

% Agents who develop symptoms
% AgentCharacteristics.Symptoms(i,1)=
% 0: no symptoms
% 1: self isolating
% 2: self isolating, going to isolate
% 3: isolating
% 4: isolating awaiting clearance test result who will leave
% 5: isolating awaiting clearance test result who will restart isolation
% AgentCharacteristics.Case(i,2)= time left in case state (if a case)
% AgentCharacteristics.Case(i,3)= time since symptom onset
AgentCharacteristics.Symptoms = zeros(parameters.PopSize, 3);

% Agents who are contact traced:
% AgentCharacteristics.Quarantine(i,1)=
% 0: not currently being traced/managed
% 1: identified, awaiting entry into quarantine
% 2: in quarantine, going to isolation
% 3: in quarantine, going to remain in quarantine
% 4: quarantining
% 5: quarantining awaiting clearance test result who will leave
% 6: quarantining awaiting clearance test result who will enter isolation
% AgentCharacteristics.Quarantine (i,2)=time left in state (if being managed)
AgentCharacteristics.Quarantine = zeros(parameters.PopSize, 2);

% Agents who are in lockdown:
% AgentCharacteristics.Lockdown(i,1)=
% 0: not currently in lockdown
% 1: identified, awaiting entry into lockdown
% 2: in lockdown, going to isolation, waiting on test result
% 3: in lockdown, going to isolation, have test result, waiting on entering isolation
% 4: in lockdown, going to remain in lockdown, waiting on test result
% 5: lockdown
% 6: lockdown, awaiting clearance test result, will leave
% AgentCharacteristics. Lockdown (i,2)=time left in lockdown state (if in lockdown)
AgentCharacteristics.Lockdown = zeros(parameters.PopSize, 2);

% AgentCharacteristics.VaccinationStatus(i,1): number of doses of Pfizer
% AgentCharacteristics.VaccinationStatus(i,2): proportional effect of VE | # of doses 
% AgentCharacteristics.VaccinationStatus(i,3): number of doses of Astra Zenica
% AgentCharacteristics.VaccinationStatus(i,4): proportional effect of VE | # of doses
% AgentCharacteristics.VaccinationStatus(i,5): Relevant when all or nothing vaccination considered.  
% With all or nothing, VE effect is either 0 or fully protective.
% = 0 if agent has no reduction in infection/sympt frac 
% = 1 if agent is fully protected (not susceptible).

% Proportional effect of VE has a range of [0,1].  E.g., 
% - if 1 dose: effective VE = proportional effect of VE * max_VE_dose 1
% - if 2 doses: effective VE = max_VE_dose 1 + proportional effect of VE * (max_VE_dose 2 - max_VE_dose 1) 
AgentCharacteristics.VaccinationStatus = zeros(parameters.PopSize,5);

% AgentCharacteristics.ReactiveVaccination stores details about timing of scheduled vaccinations and vaccine
% hesitancy
% AgentCharacteristics.ReactiveVaccination(i,1): if awaiting second dose, time until can have second dose (units days)
% AgentCharacteristics.ReactiveVaccination(i,2): if had second dose during simulation, this is 0. 
% This is set to -1 for fully vacced at start of simulation (necessary for clinical 
% pathways model to track these agents).
% ReactiveVaccination(i,3): = 1: vaccine hesitant, 0: otherwise 
AgentCharacteristics.ReactiveVaccination = zeros(parameters.PopSize,3);
AgentCharacteristics.ReactiveVaccination(:,2) = -1;

% Randomly select agents to be exposed and vaccinated in population

% Initialisation

% Choose agents to be infected, etc
AgentIndex = (1 : 1 : parameters.PopSize)';
I1 = randperm(length(AgentIndex),parameters.NE0)';
NewAgentIndex = setdiff(AgentIndex,I1);
I2 = NewAgentIndex(randperm(length(NewAgentIndex),parameters.NPI0));
NewAgentIndex = setdiff(NewAgentIndex,I2);
I3 = NewAgentIndex(randperm(length(NewAgentIndex),parameters.NI0));
NewAgentIndex = setdiff(NewAgentIndex,I3);
I4 = NewAgentIndex(randperm(length(NewAgentIndex),parameters.NPA0));
NewAgentIndex = setdiff(NewAgentIndex,I4);
I5 = NewAgentIndex(randperm(length(NewAgentIndex),parameters.NA0));
NewAgentIndex = setdiff(NewAgentIndex,I5);
I6 = NewAgentIndex(randperm(length(NewAgentIndex),parameters.NR0));

% Update agents' infection status 
AgentCharacteristics.InfectionStatus(I1,1) = 1; 
AgentCharacteristics.TestSensitivity(I1,1) = parameters.TestSensitivity(2); 
AgentCharacteristics.InfectionStatus(I2,1) = 2;  
AgentCharacteristics.TestSensitivity(I2,1) = parameters.TestSensitivity(3); 
AgentCharacteristics.InfectionStatus(I3,1) = 3; 
AgentCharacteristics.TestSensitivity(I3,1) = parameters.TestSensitivity(4); 
AgentCharacteristics.InfectionStatus(I4,1) = 4;  
AgentCharacteristics.TestSensitivity(I4,1) = parameters.TestSensitivity(5); 
AgentCharacteristics.InfectionStatus(I5,1) = 5; 
AgentCharacteristics.TestSensitivity(I5,1) = parameters.TestSensitivity(6); 
AgentCharacteristics.InfectionStatus(I6,1) = 6;  
AgentCharacteristics.TestSensitivity(I6,1) = parameters.TestSensitivity(7); 

% Store infection duration sets: latency, pre-symp infectious, sympt
% infectious
n0 = length([I1;I2;I3;I4;I5]);
n = n0*100;
lp = lognrnd(parameters.MeanDurationLatency,parameters.SDLatency,n,1);
ip = lognrnd(parameters.MeanDurationIncubation,parameters.SDIncubation,n,1);
ln = lp(lp<ip); in = ip(lp<ip);
s = lognrnd(parameters.MeanDuratioSymptoms,parameters.SDSymptoms,length(in),1);

DurationInfectionRand = [ln in-ln s];
AgentCharacteristics.InfectionStatus([I1;I2;I3;I4;I5],2:4) = DurationInfectionRand(1:n0,:);

% Choose agents to be vaccinated 
% (assume full VE effect for dose received for pre-emptive scenarios)

AgentIndex = (1 : 1 : parameters.PopSize)';
age_groups = [0 12 16 40 60 parameters.AgeDeath];

if parameters.exemplar == 3
    age_groups = [0 12 16 20:10:80 parameters.AgeDeath];
end

% For each vaccine age group
for i = 1:(length(age_groups) - 1)
    
    % Find agents in age groups
    relevant_agents = AgentIndex(AgentCharacteristics.Age >= age_groups(i) & ...
        AgentCharacteristics.Age < age_groups(i+1));
    
    % Count them up
    number_relevant_agents = length(relevant_agents);
    
    % Work out how many to vaccinate with first dose Pfizer
    number_agents_get_vaccine = round(parameters.VP10(i)*number_relevant_agents);
    
    % Select agents to be vacced and store these agents indices in V1 
    if number_agents_get_vaccine > 0
        V1 = relevant_agents(randperm(number_relevant_agents,number_agents_get_vaccine))';
    else
        V1 = [];
    end
    
    % Remove V1 from the list of agent indices in age group
    NewAgentIndex = setdiff(relevant_agents,V1);
    
    % Work out how many to vaccinate with second dose Pfizer
    number_agents_get_vaccine = round(parameters.VP20(i)*number_relevant_agents);
    
    % Select agents to be vacced and store these agents indices in V2 
    if number_agents_get_vaccine > 0
        V2 = NewAgentIndex(randperm(length(NewAgentIndex),number_agents_get_vaccine));
    else
        V2=[];
    end
    
    % Remove V2 from the list of agent indices in age group
    NewAgentIndex = setdiff(NewAgentIndex,V2);
    
    % Work out how many to vaccinate with first dose AZ
    number_agents_get_vaccine = round(parameters.VAZ10(i)*number_relevant_agents);
    
    % Select agents to be vacced and store these agents indices in V3 
    if number_agents_get_vaccine > 0
        V3 = NewAgentIndex(randperm(length(NewAgentIndex),number_agents_get_vaccine));
    else
        V3=[];
    end
    
    % Remove V2 from the list of agent indices in age group
    NewAgentIndex = setdiff(NewAgentIndex,V3);
    
    % Work out how many to vaccinate with second dose AZ
    number_agents_get_vaccine = round(parameters.VAZ20(i)*number_relevant_agents);
    if number_agents_get_vaccine > length(NewAgentIndex)
        number_agents_get_vaccine = length(NewAgentIndex);
    end
    
    % Select agents to be vacced and store these agents indices in V4 
    if number_agents_get_vaccine > 0
        V4 = NewAgentIndex(randperm(length(NewAgentIndex),number_agents_get_vaccine));
    else
        V4 = [];
    end
    
    % Store vaccination status for relevant agents
    AgentCharacteristics.VaccinationStatus(V1,1) = 1;  
    AgentCharacteristics.VaccinationStatus(V2,1) = 2; 
    AgentCharacteristics.VaccinationStatus([V1(:); V2(:)],2) = 1;   
    AgentCharacteristics.VaccinationStatus(V3,3) = 1;   
    AgentCharacteristics.VaccinationStatus(V4,3) = 2;   
    AgentCharacteristics.VaccinationStatus([V3(:); V4(:)],4) = 1; 
    
    % Reactive Vaccination.  Assume agents with one dose are scheduled to
    % get their second dose (time to dose is uniform distribution)
    if parameters.reactivevacc == 1
        numbertimestepsVE = length(parameters.VE_Effect_Over_Time_Dose1);
        % Pfizer
        TimeToSecondDose = parameters.DoseIntervalP * rand(size(V1));
        ProportionAlongTimeDependentVE = (parameters.DoseIntervalP - TimeToSecondDose)./21;
        IndVE = ceil(numbertimestepsVE*ProportionAlongTimeDependentVE);
        VE = ones(size(TimeToSecondDose));      
        VE(ProportionAlongTimeDependentVE<=1) = ...
            parameters.VE_Effect_Over_Time_Dose1(IndVE(ProportionAlongTimeDependentVE<=1));
        AgentCharacteristics.VaccinationStatus(V1,2) = VE; 
        AgentCharacteristics.ReactiveVaccination(V1,1) = TimeToSecondDose;
        % AZ
        TimeToSecondDose = parameters.DoseIntervalAZ * rand(size(V3));
        ProportionAlongTimeDependentVE = (parameters.DoseIntervalAZ - TimeToSecondDose)./21;
        IndVE = ceil(numbertimestepsVE*ProportionAlongTimeDependentVE);
        VE = ones(size(TimeToSecondDose));      
        VE(ProportionAlongTimeDependentVE<=1) = ...
            parameters.VE_Effect_Over_Time_Dose1(IndVE(ProportionAlongTimeDependentVE<=1));
        AgentCharacteristics.VaccinationStatus(V3,4) = VE; 
        AgentCharacteristics.ReactiveVaccination(V3,1) = TimeToSecondDose;
    end
    
    % If considering all or nothing model of vaccine protection, then
    % determine whether vaccinated have 0 or full effect of vaccine
    if parameters.allornothingve == 1
        % Dose 1 Pfizer (V1)
        RandomNumbers = rand(size(V1));
        IndexP1Full = V1(RandomNumbers <= parameters.Max_VE_Reduction_Infection_Pfizer_Dose1_AON);
        AgentCharacteristics.VaccinationStatus(IndexP1Full,5) = 1;
        
        % Dose 2 Pfizer (V2)
        RandomNumbers = rand(size(V2));
        IndexP2Full = V2(RandomNumbers <= parameters.Max_VE_Reduction_Infection_Pfizer_Dose2_AON);
        AgentCharacteristics.VaccinationStatus(IndexP2Full,5) = 1;
        
        % Dose 1 AZ (V3)
        RandomNumbers = rand(size(V3));
        IndexA1Full = V3(RandomNumbers <= parameters.Max_VE_Reduction_Infection_AZ_Dose1_AON);
        AgentCharacteristics.VaccinationStatus(IndexA1Full,5) = 1;
        
        % Dose 2 AZ (V4)
        RandomNumbers = rand(size(V4));
        IndexA2Full = V4(RandomNumbers <= parameters.Max_VE_Reduction_Infection_AZ_Dose2_AON);
        AgentCharacteristics.VaccinationStatus(IndexA2Full,5) = 1;
               
    end
    
end

% Find unvaccinated, non-infected agents, and assign hesitancy status
if parameters.reactivevacc == 1
    KeepIndexAllAgents = (1:parameters.PopSize)';
    NonVacced = KeepIndexAllAgents(AgentCharacteristics.VaccinationStatus(:,1)==0 & ...
    AgentCharacteristics.VaccinationStatus(:,3)==0 & ...
    AgentCharacteristics.InfectionStatus(:,1)==0);
    NonVacced(rand(size(NonVacced))>=parameters.Hesitancy) = [];
    AgentCharacteristics.ReactiveVaccination(NonVacced,3) = 1;
end


end
