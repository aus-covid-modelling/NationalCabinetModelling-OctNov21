
function [AgentCharacteristics, SummaryStatistics, parameters] = ...
    simulator(AgentCharacteristics, parameters)

    % This script simulates the agent-based model once, based on the
    % model parameters (stored in parameters) and given the inialised
    % agent characteristics (stored in AgentCharacteristics)

    % Pregenerate random numbers
    parameters = prn(parameters);

    % Initialise storage of summary statistics
    SummaryStatistics.AgeVaccIncidence = zeros(parameters.AgeDeath, 5, ceil((parameters.Ntimesteps - 1)*parameters.dt));
    SummaryStatistics.NumberInfectionVaccinationAgeStatusTime = zeros(7,5,parameters.AgeDeath, parameters.Ntimesteps);
    SummaryStatistics.InfectorInfecteePairs = [];
    SummaryStatistics.ClinicalPathwaysData = [];
    SummaryStatistics.TimeSymptomOnset = [];
    SummaryStatistics.NumberSusceptibleInDwellingsTime = [];
    SummaryStatistics.QuarantinePersonDays = [];
    SummaryStatistics.NumberInQTime = zeros(parameters.Ntimesteps,1);
    SummaryStatistics.NumberInLTime = zeros(parameters.Ntimesteps,1);
    SummaryStatistics.NumberInITime = zeros(parameters.Ntimesteps,1);
    SummaryStatistics.CaseData = [];
    SummaryStatistics.TestDate = [];
    
    % Calculate initial summary statistics
    SummaryStatistics = generate_summary_statistics(SummaryStatistics,1,parameters,AgentCharacteristics);

    % Time loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    first_case = 0;
    
    for i = 1 : parameters.Ntimesteps - 1 % for each time step
        
        timestep = i;

        CurrentDay = ceil(i*parameters.dt); %this is nth day of the simulation      
        KeepIndexAllAgents = 1 : length(AgentCharacteristics.ID); % this is the index of individuals

        % Store current agents' characteristics here
        CurrentInfectionStatus = AgentCharacteristics.InfectionStatus;
        CurrentVaccinationStatus = AgentCharacteristics.VaccinationStatus;
        CurrentHousehold = AgentCharacteristics.CurrentHousehold;
        CurrentSymptomStatus = AgentCharacteristics.Symptoms;
        CurrentTestSensitivity = AgentCharacteristics.TestSensitivity;
        CurrentContacts = AgentCharacteristics.Contacts;
        CurrentQuarantine = AgentCharacteristics.Quarantine;
        CurrentLockdown = AgentCharacteristics.Lockdown;
        CurrentIDs = AgentCharacteristics.ID;
        CurrentReactiveVaccination = AgentCharacteristics.ReactiveVaccination;
        CurrentAge = AgentCharacteristics.Age;
        
        % If no infected agents left, stop simulation
        InfectionStatesExcludingR = CurrentInfectionStatus(:,1);
        InfectionStatesExcludingR = InfectionStatesExcludingR(InfectionStatesExcludingR~=6);
        if sum(InfectionStatesExcludingR) == 0
            break 
        end
        
        % E -> PI/PA transition 
        [AgentCharacteristics,parameters] = endlatency(parameters,AgentCharacteristics,CurrentInfectionStatus,CurrentVaccinationStatus);
        
        % PI -> I, PA -> A transition
        [AgentCharacteristics,parameters,SummaryStatistics] = endincubation(parameters,AgentCharacteristics,CurrentInfectionStatus,SummaryStatistics,i+1);
        
        % I -> R, A -> R transition
        [AgentCharacteristics,parameters] = endinfection(parameters,AgentCharacteristics,CurrentInfectionStatus);

        % Transmission
        [AgentCharacteristics, parameters, SummaryStatistics] = ...
            transmission(parameters, AgentCharacteristics, CurrentInfectionStatus, ...
            CurrentSymptomStatus, CurrentQuarantine, CurrentLockdown, ...
            CurrentVaccinationStatus, CurrentReactiveVaccination, CurrentHousehold, KeepIndexAllAgents, ...
            SummaryStatistics, CurrentDay,i+1);

        % Household mobility (occurs once per day)
        if mod(i*parameters.dt,1)==0
            [parameters, AgentCharacteristics] = hhmobility(parameters,AgentCharacteristics,CurrentHousehold,CurrentSymptomStatus, CurrentQuarantine, CurrentLockdown);
        end
                
        % Update vaccine efficacy and give second doses for agents with one dose
        [AgentCharacteristics, parameters] = update_vaccine_efficacy(parameters,AgentCharacteristics,...
            CurrentVaccinationStatus,KeepIndexAllAgents,CurrentReactiveVaccination);
        
        % If response is switched on
        if parameters.response == 1
        
            % Manage agents who develop symptoms
            [AgentCharacteristics, SummaryStatistics, parameters] = manage_agents_with_symptoms(parameters, AgentCharacteristics, SummaryStatistics, ...
                CurrentSymptomStatus,CurrentTestSensitivity, CurrentContacts, CurrentInfectionStatus, KeepIndexAllAgents, CurrentIDs, i+1, CurrentHousehold, CurrentAge, CurrentVaccinationStatus);

            % Manage agents in quarantine
            [AgentCharacteristics, SummaryStatistics, parameters] = manage_quarantine(parameters, AgentCharacteristics, SummaryStatistics, ...
                CurrentSymptomStatus, CurrentTestSensitivity, CurrentContacts, CurrentInfectionStatus, CurrentQuarantine, KeepIndexAllAgents, CurrentIDs, i+1, CurrentHousehold, CurrentAge, CurrentVaccinationStatus);
           
            % If lockdown switched on
            if parameters.lockdown == 1
                
                % Enact lockdown once first case detected
                if first_case == 0
                    time_remaining = parameters.MinDurationLoc + parameters.DelayCaseIdentifiedLockdown;
                    AgentsIsolating = CurrentSymptomStatus(:,1);
                    AgentsIsolating = AgentsIsolating(AgentsIsolating==2);
                    if ~isempty(AgentsIsolating)
                        first_case = 1;  
                        AgentCharacteristics.Lockdown(:,1) = 1;
                        AgentCharacteristics.Lockdown(:,2) = parameters.DelayCaseIdentifiedLockdown;
                    end
                elseif first_case == 1
                    [AgentCharacteristics, SummaryStatistics, parameters] = ...
                            lockdown(parameters,AgentCharacteristics,SummaryStatistics,...
                            CurrentSymptomStatus, CurrentTestSensitivity,CurrentContacts, CurrentInfectionStatus,...
                            CurrentLockdown, KeepIndexAllAgents, CurrentIDs, i+1, CurrentHousehold, CurrentAge, CurrentVaccinationStatus);
                    
                    time_remaining = time_remaining - parameters.dt;
                    
                    % Turn off lockdown if duration of lockdown exceeds
                    % paramesters.MinDurationLoc
                    if time_remaining <= 0
                        parameters.lockdown = 0;
                        AgentCharacteristics.Lockdown(:,1) = 0;
                    end
                end
               
            end

            % If reactive vaccination switched on
            if parameters.reactivevacc == 1   
                % Enact program once first case detected
                if first_case == 0
                    time_to_start_rv = parameters.DelayStartReactiveVaccination;
                    AgentsIsolating = CurrentSymptomStatus(:,1);
                    AgentsIsolating = AgentsIsolating(AgentsIsolating==2);
                    if ~isempty(AgentsIsolating)
                        first_case = 1;  
                    end
                elseif first_case == 1
                    time_to_start_rv = time_to_start_rv - parameters.dt;
                    if time_to_start_rv <= 0
                        
                        % Only vaccinate once per day
                        if mod(i*parameters.dt,1)==0
 
                            [AgentCharacteristics, SummaryStatistics, parameters] = ...
                                reactive_vaccination(parameters,AgentCharacteristics,SummaryStatistics,...
                                CurrentInfectionStatus, CurrentVaccinationStatus, KeepIndexAllAgents', ...
                                CurrentReactiveVaccination);

                        end
                    end
                end
            end
        end
        
        % Update test sensitivity
        AgentCharacteristics.Symptoms(:,3) = AgentCharacteristics.Symptoms(:,3) + parameters.dt;
        
        % Aging
        [AgentCharacteristics, parameters] = aging(parameters,AgentCharacteristics);

        % Calculate summary statistics
        SummaryStatistics = generate_summary_statistics(SummaryStatistics,i+1,parameters,AgentCharacteristics);
        
    end
    
    % End time loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Collect data for clinical pathways model
    if ~isempty(SummaryStatistics.ClinicalPathwaysData)
        [a,ind] = ismember(SummaryStatistics.ClinicalPathwaysData(:,1),SummaryStatistics.TimeSymptomOnset(:,1));
        CPDreduced = SummaryStatistics.ClinicalPathwaysData(a,:);
        ind = ind(ind>0);
        TSOCD = SummaryStatistics.TimeSymptomOnset(ind,2);

        AgeGroups = calculate_age_groups(CPDreduced(:,2),parameters.TPACD);
        CPDreduced = [floor(TSOCD*parameters.dt) AgeGroups(:) CPDreduced(:,3:6)];
    else
        CPDreduced = [];
    end
    
    % Collect case data
    % Determine date of onset (TSO) and add to CaseData:
    % - asymptomatic: day of test
    % - symptomatic: day of symptom onset
    if ~isempty(SummaryStatistics.CaseData)
        
        SummaryStatistics.CaseData(:,4) = floor(SummaryStatistics.CaseData(:,4)); %notification in days
        
        % Remove any repeat entries for cases found via self iso and
        % quaranting/ld
        [~,ind] = unique(SummaryStatistics.CaseData(:,1));
        SummaryStatistics.CaseData = SummaryStatistics.CaseData(ind,:);

        IDcases = SummaryStatistics.CaseData(:,1);
        
        % symptomatic
        [logics, ind] = ismember(IDcases,SummaryStatistics.TimeSymptomOnset(:,1));
        TSOcasesS = floor(SummaryStatistics.TimeSymptomOnset(ind(ind>0),2)*parameters.dt); %TSO in timesteps
        CaseData = [SummaryStatistics.CaseData(logics,:) TSOcasesS];
        
        % asymptomatic
        IDcasesA = IDcases(~logics);
        [logics, ind] = ismember(IDcasesA,SummaryStatistics.TestDate(:,1));
        TSOcasesA = floor(SummaryStatistics.TestDate(ind(ind>0),2)); %test data in days
        CaseData = [CaseData;
        SummaryStatistics.CaseData(logics,:) TSOcasesA];
        
    else
        CaseData = [];
    end
    
    % Reduce size of Summary Statistics by only collecting daily data 
    if parameters.response == 1
    
        ntpd = 1/parameters.dt; % number of timesteps per day
        timepoints1 = 1:ntpd:timestep; % these are the relevant timesteps
        if timestep>=parameters.Ntimesteps-1
            timepoints1 = 1:ntpd:(parameters.Ntimesteps-1);
        end

        if ~isempty(timepoints1)

            SSNew.TotalInfectedEachAgeVaccGroup = SummaryStatistics.AgeVaccIncidence(:,:,1:floor(timestep*parameters.dt));
            SSNew.NumberInfectionVaccinationAgeStatusTime = SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(:,:,:,timepoints1);
            SSNew.QuarantinePersonDays = SummaryStatistics.QuarantinePersonDays;
            SSNew.NumberInQTime = SummaryStatistics.NumberInQTime(timepoints1,1);
            SSNew.NumberInLTime = SummaryStatistics.NumberInLTime(timepoints1,1);
            SSNew.NumberInITime = SummaryStatistics.NumberInITime(timepoints1,1);
            SSNew.ClinicalPathwaysData = CPDreduced;
            SSNew.CaseData = CaseData;
            SummaryStatistics = SSNew;
            
        else 
            SummaryStatistics = [];
        end
    end
    
    % Remove large arrays of random numbers from outputs
    parameters.UniformRand=[];
    parameters.ContactsNumberRand = [];
    parameters.DurationInfectionRand = [];

end

% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the probability of agents having symptoms, given age and
% vaccination status of infected agent
function  ProbabilitySymptoms = calculate_PSymptoms(Age,VacStatus,parameters)
    
    % This is where we will store effect of vaccination on symptoms (= 0 
    % if not vaccinated)
    VE = zeros(length(Age),1);
    
    % Find indices of vaccinated agents by P/AZ and Dose
    IndexP1 = find(VacStatus(:,1)==1);
    IndexP2 = find(VacStatus(:,1)==2);
    IndexA1 = find(VacStatus(:,3)==1);
    IndexA2 = find(VacStatus(:,3)==2);
    
    % If leaky vaccine model considered, find current level of protection
    % for vaccinated agents
    if parameters.allornothingve == 0
    
        VE(IndexP1) = parameters.Max_VE_Reduction_SymptomaticInfection_Pfizer_Dose1 .*  ...
            VacStatus(IndexP1,2);

        VE(IndexP2) = parameters.Max_VE_Reduction_SymptomaticInfection_Pfizer_Dose1 +  ...
            (parameters.Max_VE_Reduction_SymptomaticInfection_Pfizer_Dose2 - ...
            parameters.Max_VE_Reduction_SymptomaticInfection_Pfizer_Dose1) .*  ...
            VacStatus(IndexP2,2);

        VE(IndexA1) = parameters.Max_VE_Reduction_SymptomaticInfection_AZ_Dose1 .*  ...
            VacStatus(IndexA1,4);

        VE(IndexA2) = parameters.Max_VE_Reduction_SymptomaticInfection_AZ_Dose1 + ...
            (parameters.Max_VE_Reduction_SymptomaticInfection_AZ_Dose2 - ...
            parameters.Max_VE_Reduction_SymptomaticInfection_AZ_Dose1) .*  ...
            VacStatus(IndexA2,4);
    else
        % If all-or-nothing vaccine model considered, determine those who 
        % are fully protected, or have no protection
        VE(IndexP1) = VacStatus(IndexP1,5);
        VE(IndexP2) = VacStatus(IndexP2,5);
        VE(IndexA1) = VacStatus(IndexA1,5);
        VE(IndexA2) = VacStatus(IndexA2,5);
         
    end
     
    AgeGroupSymptoms = calculate_age_groups(Age,parameters.SymptomsACD);
    % Determine age effect on deveopment of symptoms
    ProbabilitySymptoms = parameters.ProbabilitySymptoms(AgeGroupSymptoms);
    % Modify for vaccinated
    ProbabilitySymptoms = ProbabilitySymptoms .* (1 - VE);
    
end

% Calculate the relative transmission probabilities (RelativeTP) given age 
% and vaccination status of contacts (Age, VacStatus), and relative 
% infectiousness of infected to various age groups (TPinfected)
function  RelativeTP = calculate_tp(Age,VacStatus,parameters,TPinfected)
    
    % This is where we will store effect of vaccination on transmission 
    % (= 0 if contact not vaccinated)
    VE = zeros(length(Age),1);
 
    % Determine contacts that are vaccinated
    IndexP1 = find(VacStatus(:,1)==1);
    IndexP2 = find(VacStatus(:,1)==2);
    IndexA1 = find(VacStatus(:,3)==1);
    IndexA2 = find(VacStatus(:,3)==2);
    
    % If leaky vaccine model considered, find current level of protection
    % for vaccinated agents
    if parameters.allornothingve == 0
    
        VE(IndexP1) = parameters.Max_VE_Reduction_Infection_Pfizer_Dose1 .*  ...
            VacStatus(IndexP1,2);

        VE(IndexP2) = parameters.Max_VE_Reduction_Infection_Pfizer_Dose1 +  ...
            (parameters.Max_VE_Reduction_Infection_Pfizer_Dose2 - ...
            parameters.Max_VE_Reduction_Infection_Pfizer_Dose1) .*  ...
            VacStatus(IndexP2,2);

        VE(IndexA1) = parameters.Max_VE_Reduction_Infection_AZ_Dose1 .*  ...
            VacStatus(IndexA1,4);

        VE(IndexA2) = parameters.Max_VE_Reduction_Infection_AZ_Dose1 + ...
            (parameters.Max_VE_Reduction_Infection_AZ_Dose2 - ...
            parameters.Max_VE_Reduction_Infection_AZ_Dose1) .*  ...
            VacStatus(IndexA2,4);
    else
        % If all-or-nothing vaccine model considered, determine those who 
        % are fully protected, or have no protection
        VE(IndexP1) = VacStatus(IndexP1,5);
        VE(IndexP2) = VacStatus(IndexP2,5);
        VE(IndexA1) = VacStatus(IndexA1,5);
        VE(IndexA2) = VacStatus(IndexA2,5);
        
    end
    
    AgeGroupContacts = calculate_age_groups(Age,parameters.TPACD);
    % Determine age effect on susceptibility
    RelativeTP = TPinfected(1,AgeGroupContacts);
    % Determine overall effect on relative transmission probabilities, 
    % incorporating age and vacc status of contacts, and infectiousness of 
    % infected 
    RelativeTP = RelativeTP(:) .* (1 - VE(:));

end

% Update agents due to transition E -> PI/PA  
function [AgentCharacteristics,parameters] = endlatency(parameters,AgentCharacteristics,...
    CurrentInfectionStatus,CurrentVaccinationStatus)
    
    % Find index of agents that are in latent phase 
    IndexAgents = find(CurrentInfectionStatus(:,1) == 1);
    IndexAgentsCopy = IndexAgents;
    
    % For these agents, find those that have timer>0, and remove from list
    IndexAgents(CurrentInfectionStatus(IndexAgents,2) > 0) = [];
    
    % Keep list of other agents for updating latency timer:
    IndexAgentsCopy(CurrentInfectionStatus(IndexAgents,2) <= 0) = [];
    
    % Calculate the probability of transitioning agents having symptoms
    ProbabilitySymptoms = calculate_PSymptoms(AgentCharacteristics.Age(IndexAgents),...
        CurrentVaccinationStatus(IndexAgents,:),parameters);
    
    % Determine agents that have symptoms and update infection and test 
    % sensitivity status
    RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexAgents)-1);
    IndexAgentsPA = IndexAgents(RandomNumbers > ProbabilitySymptoms);
    IndexAgentsPI = IndexAgents(RandomNumbers <= ProbabilitySymptoms);
    
    % Update counter for pregenerated random numbers:
    parameters.countUR = parameters.countUR + length(IndexAgents);
    if parameters.countUR > 0.9*10^5
        parameters = updatecounter(parameters, "uniform");
    end

    AgentCharacteristics.InfectionStatus(IndexAgentsPA,1) = 4;
    AgentCharacteristics.TestSensitivity(IndexAgentsPA,1) = parameters.TestSensitivity(5);
    AgentCharacteristics.InfectionStatus(IndexAgentsPI,1) = 2;
    AgentCharacteristics.TestSensitivity(IndexAgentsPI,1) = parameters.TestSensitivity(3);
    
    % Update latency timers for transitioning agents
    AgentCharacteristics.InfectionStatus(IndexAgentsPA,2) = 0;
    AgentCharacteristics.InfectionStatus(IndexAgentsPI,2) = 0;
    
    % Update latency timers for non-transitioning agents
    AgentCharacteristics.InfectionStatus(IndexAgentsCopy,2) = ...
        AgentCharacteristics.InfectionStatus(IndexAgentsCopy,2) - parameters.dt;
    
end


% Update agents due to transition PI/PA -> I/A  
function [AgentCharacteristics,parameters,SummaryStatistics] = endincubation(parameters,...
    AgentCharacteristics,CurrentInfectionStatus,SummaryStatistics,timestep)
      
    % Find index of agents that are presymptomatic  
    IndexAgentsPI = find(CurrentInfectionStatus(:,1) == 2);
    IndexAgentsPA = find(CurrentInfectionStatus(:,1) == 4);
    IndexAgentsPICopy = IndexAgentsPI;
    IndexAgentsPACopy = IndexAgentsPA;
    
    % For these agents, find those that have timer>0, and remove from list
    IndexAgentsPI(CurrentInfectionStatus(IndexAgentsPI,3) > 0) = [];
    IndexAgentsPA(CurrentInfectionStatus(IndexAgentsPA,3) > 0) = [];
    
    % Keep list of other agents for updating presymptomatic timer:
    IndexAgentsPICopy(CurrentInfectionStatus(IndexAgentsPI,3) <= 0) = [];
    IndexAgentsPACopy(CurrentInfectionStatus(IndexAgentsPA,3) <= 0) = [];
    
    % Update infection and test sensitivity status for transitioning agents 
    AgentCharacteristics.InfectionStatus(IndexAgentsPI,1) = 3;
    AgentCharacteristics.TestSensitivity(IndexAgentsPI,1) = parameters.TestSensitivity(4);
    AgentCharacteristics.InfectionStatus(IndexAgentsPA,1) = 5;
    AgentCharacteristics.TestSensitivity(IndexAgentsPA,1) = parameters.TestSensitivity(6);
    
    % Store time of symptom onset for symptomatic agents
    SummaryStatistics.TimeSymptomOnset = [SummaryStatistics.TimeSymptomOnset;
        AgentCharacteristics.ID(IndexAgentsPI) timestep * ones(size(IndexAgentsPI))];
    
    % Update isolation pathway details
    AgentCharacteristics.Symptoms(IndexAgentsPI,1) = 1;
    AgentCharacteristics.Symptoms(IndexAgentsPI,2) = parameters.DelayOnsetSymptTestResult;
    AgentCharacteristics.Symptoms(IndexAgentsPI,3) = 0;
    AgentCharacteristics.Symptoms(IndexAgentsPA,3) = 0;
    
    % Update presymptomatic timers for transitioning agents
    AgentCharacteristics.InfectionStatus([IndexAgentsPI; IndexAgentsPA],3) = 0;
    
    % Update presymptomatic timers for non-transitioning agents
    AgentCharacteristics.InfectionStatus([IndexAgentsPICopy; IndexAgentsPACopy],3) = ...
        AgentCharacteristics.InfectionStatus([IndexAgentsPICopy; IndexAgentsPACopy],3) - parameters.dt;         
    
end

% Update agents due to transition I/A -> R
function [AgentCharacteristics,parameters] = endinfection(parameters,...
    AgentCharacteristics,CurrentInfectionStatus)
      
    % Find index of agents that are symptomatic  
    IndexAgentsI = find(CurrentInfectionStatus(:,1) == 3);
    IndexAgentsA = find(CurrentInfectionStatus(:,1) == 5);
    IndexAgentsICopy = IndexAgentsI;
    IndexAgentsACopy = IndexAgentsA;
    
    % For these agents, find those that have timer>0, and remove from list
    IndexAgentsI(CurrentInfectionStatus(IndexAgentsI,4) > 0) = [];
    IndexAgentsA(CurrentInfectionStatus(IndexAgentsA,4) > 0) = [];
    
    % Keep list of other agents for updating symptomatic timer:
    IndexAgentsICopy(CurrentInfectionStatus(IndexAgentsI,4) <= 0) = [];
    IndexAgentsACopy(CurrentInfectionStatus(IndexAgentsA,4) <= 0) = [];
    
    % Update infection status for transitioning agents 
    AgentCharacteristics.InfectionStatus([IndexAgentsI; IndexAgentsA],1) = 6;
    AgentCharacteristics.TestSensitivity([IndexAgentsI; IndexAgentsA],1) = parameters.TestSensitivity(7);
    
    % Update symptomatic timers for transitioning agents
    AgentCharacteristics.InfectionStatus([IndexAgentsI; IndexAgentsA],4) = 0;
    
    % Update symptomatic timers for non-transitioning agents
    AgentCharacteristics.InfectionStatus([IndexAgentsICopy; IndexAgentsACopy],4) = ...
        AgentCharacteristics.InfectionStatus([IndexAgentsICopy; IndexAgentsACopy],4) - parameters.dt;
    
end

% Update vaccine efficacy and give second doses for agents with one dose
function [AgentCharacteristics, parameters] = update_vaccine_efficacy(parameters,...
    AgentCharacteristics,CurrentVaccinationStatus, KeepIndexAllAgents, CurrentReactiveVaccination)
    
    % First, update vaccine efficacy
    
    %  Find row indices of agents that are vaccinated and have not reached
    %  full effect of vaccine dose
    IndexAgents = KeepIndexAllAgents'; % all agents
    IndexP1 = IndexAgents(CurrentVaccinationStatus(:,1)==1 & CurrentVaccinationStatus(:,2)~=1);
    IndexP2 = IndexAgents(CurrentVaccinationStatus(:,1)==2 & CurrentVaccinationStatus(:,2)~=1);
    IndexA1 = IndexAgents(CurrentVaccinationStatus(:,3)==1 & CurrentVaccinationStatus(:,4)~=1);
    IndexA2 = IndexAgents(CurrentVaccinationStatus(:,3)==2 & CurrentVaccinationStatus(:,4)~=1);
    
    % Get current VEs
    VEDoseP1 = CurrentVaccinationStatus(IndexP1,2);
    VEDoseA1 = CurrentVaccinationStatus(IndexA1,4);
    VEDoseP2 = CurrentVaccinationStatus(IndexP2,2);
    VEDoseA2 = CurrentVaccinationStatus(IndexA2,4);
    
    % Update VEs for these agents
    [~, ind1p]= ismember(VEDoseP1,parameters.VE_Effect_Over_Time_Dose1);
    [~, ind1a]= ismember(VEDoseA1,parameters.VE_Effect_Over_Time_Dose1);
    [~, ind2p]= ismember(VEDoseP2,parameters.VE_Effect_Over_Time_Dose2);
    [~, ind2a]= ismember(VEDoseA2,parameters.VE_Effect_Over_Time_Dose2);
    
    VEDoseP1 = parameters.VE_Effect_Over_Time_Dose1(ind1p+1);
    VEDoseA1 = parameters.VE_Effect_Over_Time_Dose1(ind1a+1);
    VEDoseP2 = parameters.VE_Effect_Over_Time_Dose2(ind2p+1);
    VEDoseA2 = parameters.VE_Effect_Over_Time_Dose2(ind2a+1);
    
    AgentCharacteristics.VaccinationStatus(IndexP1,2) = VEDoseP1;
    AgentCharacteristics.VaccinationStatus(IndexA1,4) = VEDoseA1;
    AgentCharacteristics.VaccinationStatus(IndexP2,2) = VEDoseP2;
    AgentCharacteristics.VaccinationStatus(IndexA2,4) = VEDoseA2; 
    
    % If all-or-nothing vaccine protection model considered, instead of 
    % having proportion of full vaccine efficacy, there is instead a chance 
    % that in the the first fews following vaccination agents will
    % develop full vaccine protection.  This is determined below.
    if parameters.allornothingve == 1 

        % For agents who are not currently assigned full effect after
        % dose 1 or 2, check if they develop full protection this timestep
        
        % For dose 1, only applies to agents who are beyond the first week 
        % of receiving the dose
        IndexP1 = IndexAgents(CurrentVaccinationStatus(:,1)==1 & CurrentVaccinationStatus(:,2)~=1 & CurrentVaccinationStatus(:,2)>0.01);
        IndexA1 = IndexAgents(CurrentVaccinationStatus(:,3)==1 & CurrentVaccinationStatus(:,4)~=1 & CurrentVaccinationStatus(:,4)>0.01);
        
        % Pfizer dose 1. 
        AoN = CurrentVaccinationStatus(IndexP1,5);
        IndexP1 = IndexP1(AoN==0);
        RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexP1)-1);
        IndexP1Full = IndexP1(RandomNumbers <= parameters.ProbabilityPerTimeStep0to1VE_Pfizer_Dose1);
        % Update counter for pregenerated random numbers:
        parameters.countUR = parameters.countUR + length(IndexP1);
        if parameters.countUR > 0.9*10^5
            parameters = updatecounter(parameters, "uniform");
        end
        AgentCharacteristics.VaccinationStatus(IndexP1Full,5) = 1;
    
        % AZ dose 1
        AoN = CurrentVaccinationStatus(IndexA1,5);
        IndexA1 = IndexA1(AoN==0);
        RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexA1)-1);
        IndexA1Full = IndexA1(RandomNumbers <= parameters.ProbabilityPerTimeStep0to1VE_AZ_Dose1);   
        % Update counter for pregenerated random numbers:
        parameters.countUR = parameters.countUR + length(IndexA1);
        if parameters.countUR > 0.9*10^5
            parameters = updatecounter(parameters, "uniform");
        end
        AgentCharacteristics.VaccinationStatus(IndexA1Full,5) = 1;
        
        % Pfizer dose 2.  
        AoN = CurrentVaccinationStatus(IndexP2,5);
        IndexP2 = IndexP2(AoN==0);       
        RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexP2)-1);
        IndexP2Full = IndexP2(RandomNumbers <= parameters.ProbabilityPerTimeStep0to1VE_Pfizer_Dose2);  
        % Update counter for pregenerated random numbers:
        parameters.countUR = parameters.countUR + length(IndexP2);
        if parameters.countUR > 0.9*10^5
            parameters = updatecounter(parameters, "uniform");
        end
        AgentCharacteristics.VaccinationStatus(IndexP2Full,5) = 1;

        % AZ dose 2
        AoN = CurrentVaccinationStatus(IndexA2,5);
        IndexA2 = IndexA2(AoN==0); 
        RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexA2)-1);
        IndexA2Full = IndexA2(RandomNumbers <= parameters.ProbabilityPerTimeStep0to1VE_AZ_Dose2);   
        % Update counter for pregenerated random numbers:
        parameters.countUR = parameters.countUR + length(IndexA2);
        if parameters.countUR > 0.9*10^5
            parameters = updatecounter(parameters, "uniform");
        end
        AgentCharacteristics.VaccinationStatus(IndexA2Full,5) = 1;

    end
    
    % If reactive vaccination is switched on, then it is possible that
    % there will be single-dosed agents who will obtain second dose during 
    % simulation.  We update their vaccination status here:
    if parameters.response == 1
        if parameters.reactivevacc == 1

            % Find all agents with one dose who get second dose this timestep
            % Pfizer
            IndexP1to2 = IndexAgents(CurrentVaccinationStatus(:,1)==1 & CurrentReactiveVaccination(:,1)<=0);
            % AstraZeneca
            IndexA1to2 = IndexAgents(CurrentVaccinationStatus(:,3)==1 & CurrentReactiveVaccination(:,1)<=0);
            % Update their status
            % Pfizer
            AgentCharacteristics.VaccinationStatus(IndexP1to2,1) = 2;
            AgentCharacteristics.VaccinationStatus(IndexP1to2,2) = 0;
            AgentCharacteristics.ReactiveVaccination(IndexP1to2,2) = 0;
            % AstraZeneca
            AgentCharacteristics.VaccinationStatus(IndexA1to2,3) = 2;
            AgentCharacteristics.VaccinationStatus(IndexA1to2,4) = 0;
            AgentCharacteristics.ReactiveVaccination(IndexA1to2,2) = 0;

            % Otherwise, update timers for second dose schedule
            IndexP1 = IndexAgents(CurrentVaccinationStatus(:,1)==1 & CurrentReactiveVaccination(:,1)>0);
            IndexA1 = IndexAgents(CurrentVaccinationStatus(:,3)==1 & CurrentReactiveVaccination(:,1)>0);
            AgentCharacteristics.ReactiveVaccination(IndexP1,1) = AgentCharacteristics.ReactiveVaccination(IndexP1,1) - parameters.dt;
            AgentCharacteristics.ReactiveVaccination(IndexA1,1) = AgentCharacteristics.ReactiveVaccination(IndexA1,1) - parameters.dt;

        end
    end

end

% Update agents due to transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AgentCharacteristics, parameters, SummaryStatistics] = ...
        transmission(parameters, AgentCharacteristics, CurrentInfectionStatus, ...
        CurrentSymptomStatus, CurrentQuarantine, CurrentLockdown, ...
        CurrentVaccinationStatus, CurrentReactiveVaccination, CurrentHousehold, ...
        KeepIndexAllAgents, SummaryStatistics, CurrentDay, timestep)

    % Find row indices of agents that are infectious and susceptible. 
    IndexSymptomatic = KeepIndexAllAgents'; % all agents
    IndexAsymptomatic = KeepIndexAllAgents'; % all agents
    IndexSusceptible = KeepIndexAllAgents'; % all agents
    
    IndexSymptomatic(CurrentInfectionStatus(:,1) < 2 | CurrentInfectionStatus(:,1) > 3) = []; % this includes PI and I
    IndexAsymptomatic(CurrentInfectionStatus(:,1) < 4 | CurrentInfectionStatus(:,1) > 5) = []; % this includes PA and A
    IndexSusceptible(CurrentInfectionStatus(:,1) ~=0) = [];
    
    % Initialise storage arrays
    InfectedAgentsThisTimeStep = []; % row index of infected agents
    InfectorInfecteePairsThisTimeStep = []; % infector infectee pair details
    ClinicalPathwaysData = []; % data for clinical pathways model
    
    % Create arrays of row indices and contact rate age groups for infected
    % agents this timestep
    AgeGroupAllAgents = calculate_age_groups(AgentCharacteristics.Age,parameters.AgeClassDividersContacts);
    RowIndexInfected = [IndexSymptomatic; IndexAsymptomatic];
    AgeGroupInfected = AgeGroupAllAgents(RowIndexInfected);
    
    % Create arrays of age groups relevant for transmission probabilities
    % and VE onward transmission correction factor
    AgeGroupAllAgentsTP = calculate_age_groups(AgentCharacteristics.Age,parameters.TPACD);
    AgeGroupInfectedTP = AgeGroupAllAgentsTP(RowIndexInfected);
    
    % Calculate relative onwards transmission for vaccinated agents
    VaccineStatusInfected = CurrentVaccinationStatus(RowIndexInfected,:);
    % This is 1 for unvacced agents
    OnwardTVacc = ones(size(RowIndexInfected));
    % For vacced agents, this depends on vaccine type and number of doses
    OnwardTVacc(VaccineStatusInfected(:,1)==1) = ...
        parameters.OnwardsTCorrection(AgeGroupInfectedTP(VaccineStatusInfected(:,1)==1),1); % Pf dose 1
    OnwardTVacc(VaccineStatusInfected(:,1)==2) = ...
        parameters.OnwardsTCorrection(AgeGroupInfectedTP(VaccineStatusInfected(:,1)==2),2); % Pf dose 2
    OnwardTVacc(VaccineStatusInfected(:,3)==1) = ...
        parameters.OnwardsTCorrection(AgeGroupInfectedTP(VaccineStatusInfected(:,3)==1),3); % AZ dose 1
    OnwardTVacc(VaccineStatusInfected(:,3)==2) = ...
        parameters.OnwardsTCorrection(AgeGroupInfectedTP(VaccineStatusInfected(:,3)==2),4); % AZ dose 2
    
    % If there are infected agents this timestep then
    if ~isempty(RowIndexInfected)
        
        % Determine which agents are self-isolating, isolating, in
        % quarantine and/or in lockdown, and store row indices in these
        % arrays
        SelfIsoAgents = KeepIndexAllAgents(CurrentSymptomStatus(:,1)==1 | CurrentSymptomStatus(:,1)==2);
        IsoAgents = KeepIndexAllAgents(CurrentSymptomStatus(:,1)>=3);
        QAgents = KeepIndexAllAgents(CurrentQuarantine(:,1)>=2);
        LDAgents = KeepIndexAllAgents(CurrentLockdown(:,1)>=2);
        
        % Use these to determine relative transmission probabilities in
        % household and in community, based on Q/I/LD status 
        % (order here is important)
        ProbabilityContactableInHH = ones(size(KeepIndexAllAgents));
        ProbabilityContactableInHH(LDAgents) = 1 - parameters.EffectivenessLocH;
        ProbabilityContactableInHH(QAgents) = 1 - parameters.EffectivenessQuaH; %overwrites LD status
        ProbabilityContactableInHH(SelfIsoAgents) = 1 - parameters.EffectivenessSelfIsoH; %overwrites Q/LD status
        ProbabilityContactableInHH(IsoAgents) = 1 - parameters.EffectivenessIsoH; %overwrites Q/SI/LD status 
        ProbabilityContactableInHH=ProbabilityContactableInHH(:);

        ProbabilityContactableInComm = ones(size(KeepIndexAllAgents));
        ProbabilityContactableInComm(LDAgents) = (1 - parameters.EffectivenessLocC);
        ProbabilityContactableInComm(QAgents) = 1 - parameters.EffectivenessQuaC;
        ProbabilityContactableInComm(SelfIsoAgents) = 1 - parameters.EffectivenessSelfIsoC;
        ProbabilityContactableInComm(IsoAgents) = 1 - parameters.EffectivenessIsoC;
        ProbabilityContactableInComm=ProbabilityContactableInComm(:);

        % InfectionProbability stores the probability that each agent
        % will contract an infection following a contact. Only susceptible 
        % agents have non-zero probability of contracting an infection
        InfectionProbability = zeros(size(KeepIndexAllAgents'));
        InfectionProbability(IndexSusceptible) = 1;
                
        % Calculate relative onward transmission probability from 
        % infectious agents, based on their symptoms (PI/I vs PA/A)                    
        RelativeOnwardTransmission = ones(size(InfectionProbability));
        RelativeOnwardTransmission(length(IndexSymptomatic)+1:end) = ...
            RelativeOnwardTransmission(length(IndexSymptomatic)+1:end) * ...
            (1 - parameters.ReductionOnwardTAsymptomatic);
        
        % Calculate relative onward transmission from infectious agents, 
        % based on age and whether contact is in household or in community
        RelativeTPInfectedsHH = parameters.RelativeTransmissibilityHH(AgeGroupInfectedTP,:);
        RelativeTPInfectedsComm = parameters.RelativeTransmissibilityCom(AgeGroupInfectedTP,:);
        
        % For each agent that is infectious
        for j = 1 : length(RowIndexInfected)
            
            % This is their relative probability of transmission in either 
            % Community or household due to their Q/I/LD status
            PContactableC = squeeze(ProbabilityContactableInComm(RowIndexInfected(j),1));
            PContactableH = squeeze(ProbabilityContactableInHH(RowIndexInfected(j),1));
            
            % This is the household ID of the infectious agent's current
            % household
            ACH = CurrentHousehold(RowIndexInfected(j));

            % Here are the indices of other agents in the same household:
            HouseholdContacts = KeepIndexAllAgents(CurrentHousehold==ACH);

            % Here are the indices of agents outside hh:
            NonHouseholdMembers = KeepIndexAllAgents(CurrentHousehold~=ACH);
            
            % Remove infected agent from HH contact list
            HouseholdContacts(HouseholdContacts==RowIndexInfected(j))=[];
            
            % Ensure only count household contacts once per day
            if mod((timestep - 1)*parameters.dt,1)>0
                HouseholdContacts = [];
            end
            
            % Update household contacts based on Q/I/LD status of
            % infectious agent and household contacts
            RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(HouseholdContacts)-1,1);
            parameters.countUR = parameters.countUR + length(HouseholdContacts);
            if parameters.countUR > 0.9*10^5
                parameters = updatecounter(parameters, "uniform");
            end
            HouseholdContacts(RandomNumbers > PContactableH * ProbabilityContactableInHH(HouseholdContacts(:))) = [];
            
            % Now we look at non-household contacts.  Here are the age
            % groups of agents outside of infectious agents household
            AgeGroupsNonHouseholdMembers = AgeGroupAllAgents(NonHouseholdMembers);

            % Determine number of contacts with these other agents in
            % each age group
            X = parameters.ContactsNumberRand(AgeGroupInfected(j),:,parameters.countCNR);

            % Update counter for pregenerated random numbers:
            parameters.countCNR = parameters.countCNR + 1;
            if parameters.countCNR > 0.9*10^5
                parameters = updatecounter(parameters, "contacts");
            end

            % This is where we will store row indices of contacts:
            IndexOfContacts = HouseholdContacts;
             
            % This is the unique agent ID of the infectious agent (need 
            % this for storing infector infectee pair data)
            AIDj = AgentCharacteristics.ID(RowIndexInfected(j));
            
            % Determine which agents are contacted in the community, they
            % will be stored in this array
            IndCommContacts = [];
            
            % For each age group relevant for contacts
            for m = 1 : length(X)

                % If the agent makes contact with other agents, find the
                % indices of these agents, add to list of contacts
                if X(m) > 0

                    % Find all contactable agents in this age group.
                    IndexContactableAgents = NonHouseholdMembers(AgeGroupsNonHouseholdMembers==m);
                    
                    % Sample with replacement X(m) contacts
                    % from this list using pregenerated random
                    % numbers            
                    IndexOfContactsTemp = ceil(length(IndexContactableAgents) * ...
                        parameters.UniformRand(parameters.countUR:parameters.countUR+X(m)-1,1));
                    
                    % Transform indices to original indices and include in
                    % the list of indices of contacts
                    try
                      tempca = IndexContactableAgents(IndexOfContactsTemp);
                    catch
                      tempca = [];
                    end
                    
                    IndCommContacts = [IndCommContacts(:); tempca(:)];
                    
                    % Update counter for pregenerated random numbers
                    parameters.countUR = parameters.countUR + X(m);
                    if parameters.countUR > 0.9*10^5
                        parameters = updatecounter(parameters, "uniform");
                    end

                end
                
            end
            
            % Update community contacts based on Q/I/LD status of
            % infectious agent and community contacts
            RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndCommContacts)-1,1);
            parameters.countUR = parameters.countUR + length(IndCommContacts);
            if parameters.countUR > 0.9*10^5
                parameters = updatecounter(parameters, "uniform");
            end
            
            IndCommContacts(RandomNumbers > PContactableC * ProbabilityContactableInComm(IndCommContacts(:))) = [];
            
            % Calculate relative probability of transmission, given age of
            % infected, and age and vaccine status of contacts 
            RelativeTPComm = ...
                calculate_tp(AgentCharacteristics.Age(IndCommContacts),...
                    CurrentVaccinationStatus(IndCommContacts,:),parameters, RelativeTPInfectedsComm(j,:));
            RelativeTPHH = ...
                calculate_tp(AgentCharacteristics.Age(IndexOfContacts),...
                    CurrentVaccinationStatus(IndexOfContacts,:),parameters, RelativeTPInfectedsHH(j,:));
         
            RelativeTP = [RelativeTPHH(:); RelativeTPComm(:)];
            
            IndexOfContacts = [IndexOfContacts IndCommContacts'];
            
            % Now we calculate actual transmission probabilities, IP, for  
            % each contact.  This is:
            % IP = (for each contact, 1/0 if susc/not susc) * 
            % (for each contact, rel susc * (1-VEinf)) *
            % beta * (rel onward trans based on symptoms of infected) * 
            % (onwards correction factor for infected)
            IP = InfectionProbability(IndexOfContacts) .* RelativeTP .* ...
                parameters.Ptransmission * RelativeOnwardTransmission(j) * ...
                OnwardTVacc(j);
            
            HHContacts = zeros(size(IndexOfContacts)); % this is for storing whether contacts occur within hh
            HHContacts(1:length(HouseholdContacts))=1; % = 1 if hh contact, 0 otherwise
            
            % Add contacts to infectious agent's contact list (for tracing)
            AgentCharacteristics.Contacts{RowIndexInfected(j),1} = ...
                [AgentCharacteristics.Contacts{RowIndexInfected(j),1};
                AgentCharacteristics.ID(IndexOfContacts(:)) HHContacts(:) timestep*ones(length(IndexOfContacts),1)];
                
            % For each contact, determine whether transmission occurs.
            % Transmission Rule: More than one transmission event can occur 
            % in the time step, but each susceptible hosts may only acquire 
            % up to one infection per time step

            if ~isempty(IndexOfContacts)

                SusceptibilityHHContacts = InfectionProbability(HouseholdContacts);
                NumbersuceptibleHHcontacts = sum(SusceptibilityHHContacts>0);

                % Determine which transmissions are successful.
                % NewInfections = 1 if successful, 0 unsuccessful
                % for each contact event
                NewInfections = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexOfContacts)-1,1)<IP;
              
                % Update counter for pregenerated random numbers:
                parameters.countUR = parameters.countUR + length(IndexOfContacts);
                if parameters.countUR > 0.9*10^5
                    parameters = updatecounter(parameters, "uniform");
                end

                if any(NewInfections)

                    % Need to check transmission rules.  Does anyone
                    % acquire more than one infection?
                    % Find new infections:
                    IndexOfContacts = IndexOfContacts(NewInfections == 1);
                    HHContacts = HHContacts(NewInfections == 1);

                    % If one contact is infected more than once, remove
                    % extra successful transmissions.  Do this by
                    % finding duplicate contact indices and
                    % removing extras from IndexOfContacts
                    % SLOW:
                    % [IndexOfContacts, ~] = unique(IndexOfContacts);
                    % FAST:
                    [IndexOfContacts_temp, idx] = sort(IndexOfContacts);
                    IndexOfContacts = IndexOfContacts_temp([true;diff(IndexOfContacts_temp(:))>0]);

                    HHContacts = HHContacts(idx);
                    HHContacts = HHContacts([true;diff(IndexOfContacts_temp(:))>0]);

                    % Now update infection status:
                    AgentCharacteristics.InfectionStatus(IndexOfContacts,1) = 1;
                    AgentCharacteristics.TestSensitivity(IndexOfContacts,1) = parameters.TestSensitivity(2);
                    
                    % Assign durations for infection:
                    AgentCharacteristics.InfectionStatus(IndexOfContacts,2:4) = ...
                        parameters.DurationInfectionRand(parameters.countDIR : parameters.countDIR + length(IndexOfContacts) - 1,:);

                    % Update counter for pregenerated random numbers:
                    parameters.countDIR = parameters.countDIR + length(IndexOfContacts);
                    if parameters.countDIR > 0.9*size(parameters.DurationInfectionRand,1)
                        parameters = updatecounter(parameters, "durationinfection");
                    end

                    % Make infected susceptibles not susceptible for next
                    % infected agent of interest:
                    InfectionProbability(IndexOfContacts) = 0;

                    InfectedAgentsThisTimeStep = [InfectedAgentsThisTimeStep;IndexOfContacts'];
                    
                    % Find number of susceptible agents in households of
                    % all infected contacts (this counts agents that are now just infected)
                    HHallcontacts = AgentCharacteristics.CurrentHousehold(IndexOfContacts(:));
                    NumbersuceptibleInHHAllcontacts = SummaryStatistics.NumberSusceptibleInDwellingsTime(HHallcontacts,timestep-1);

                    % Store infectee/infector pairs 
                    % [ID_infector ID_infectee Age_infectee HH_contact HHID_contact NumSusceptible_contact timestep]
                    % HH_contact = 1 if transmission is due to a household contact
                    InfectorInfecteePairsThisTimeStep = [InfectorInfecteePairsThisTimeStep;
                        AIDj * ones(length(IndexOfContacts),1),... 
                        AgentCharacteristics.ID(IndexOfContacts(:)),... 
                        AgentCharacteristics.Age(IndexOfContacts(:)),...
                        HHContacts(:),...
                        AgentCharacteristics.CurrentHousehold(IndexOfContacts(:)),...
                        timestep * ones(length(IndexOfContacts),1),...
                        NumbersuceptibleHHcontacts * ones(length(IndexOfContacts),1),...
                        NumbersuceptibleInHHAllcontacts(:)];
                    
                    
                    % Modify Vacc Status recorded if vaccine received too
                    % recently (for clinical pathways model)
                    CVSContact = CurrentVaccinationStatus(IndexOfContacts(:),:);
                    CRVSContact = CurrentReactiveVaccination(IndexOfContacts(:),2);
                    ModLabelsP = zeros(size(CRVSContact));
                    ModLabelsA = ModLabelsP;
                    CVSContactP = CVSContact(:,1);
                    CVSContactA = CVSContact(:,3);
                    
                    CVSContactP(CVSContactP==1 & CVSContact(:,2)<parameters.MinVEForSeverityDose1) = 0;
                    CVSContactA(CVSContactA==1 & CVSContact(:,4)<parameters.MinVEForSeverityDose1) = 0;
                    
                    % if too recently vacced, remove status
                    CVSContactP(CVSContactP==2 & CVSContact(:,2)<parameters.MinVEForSeverityDose2) = 0;
                    CVSContactA(CVSContactA==2 & CVSContact(:,4)<parameters.MinVEForSeverityDose2) = 0;
                    
                    % if 2nd vacced during sim and inf late enough, label
                    % AZ ModD1, not P ModD2
                    ModLabelsA(CVSContactA==2 & CVSContact(:,4)>=parameters.MinVEForSeverityDose2 & CRVSContact==0) = 2;
                    ModLabelsP(CVSContactP==2 & CVSContact(:,2)>=parameters.MinVEForSeverityDose2 & CRVSContact==0) = 2;
                    
                    CVSContactP(CVSContactP==2 & CVSContact(:,2)>=parameters.MinVEForSeverityDose2 & CRVSContact==0) = 0;
                    CVSContactA(CVSContactA==2 & CVSContact(:,4)>=parameters.MinVEForSeverityDose2 & CRVSContact==0) = 0;
                    
                    ClinicalPathwaysData = [ClinicalPathwaysData;
                        AgentCharacteristics.ID(IndexOfContacts(:)),... 
                        AgentCharacteristics.Age(IndexOfContacts(:)),...
                        CVSContactP(:) CVSContactA(:) ModLabelsA(:) ModLabelsP(:)];
                
                
                end
            end
        end
    end
    
    if length(InfectedAgentsThisTimeStep) ~= length(unique(InfectedAgentsThisTimeStep))
        msg = 'More than 1 infection of agent during timestep';
        error(msg)
    end
    
    % Update Summary statistics
    
    AgesInfecteds = ceil(AgentCharacteristics.Age(InfectedAgentsThisTimeStep));    
    PVaccStatusInfecteds = AgentCharacteristics.VaccinationStatus(InfectedAgentsThisTimeStep,1);
    AVaccStatusInfecteds = AgentCharacteristics.VaccinationStatus(InfectedAgentsThisTimeStep,3);
    
    AI = AgesInfecteds(PVaccStatusInfecteds==1);
    NumInfectedEachAgeGroup = histcounts(AI,0:(parameters.AgeDeath));
    SummaryStatistics.AgeVaccIncidence(:,1, CurrentDay) = ...
        SummaryStatistics.AgeVaccIncidence(:,1, CurrentDay) + NumInfectedEachAgeGroup(:);
    
    AI = AgesInfecteds(PVaccStatusInfecteds==2);
    NumInfectedEachAgeGroup = histcounts(AI,0:(parameters.AgeDeath));
    SummaryStatistics.AgeVaccIncidence(:,2, CurrentDay) = ...
        SummaryStatistics.AgeVaccIncidence(:,2, CurrentDay) + NumInfectedEachAgeGroup(:);
    
    AI = AgesInfecteds(AVaccStatusInfecteds==1);
    NumInfectedEachAgeGroup = histcounts(AI,0:(parameters.AgeDeath));
    SummaryStatistics.AgeVaccIncidence(:,3, CurrentDay) = ...
        SummaryStatistics.AgeVaccIncidence(:,3, CurrentDay) + NumInfectedEachAgeGroup(:);
    
    AI = AgesInfecteds(AVaccStatusInfecteds==2);
    NumInfectedEachAgeGroup = histcounts(AI,0:(parameters.AgeDeath));
    SummaryStatistics.AgeVaccIncidence(:,4, CurrentDay) = ...
        SummaryStatistics.AgeVaccIncidence(:,4, CurrentDay) + NumInfectedEachAgeGroup(:);
    
    AI = AgesInfecteds(AVaccStatusInfecteds==0 & PVaccStatusInfecteds==0);
    NumInfectedEachAgeGroup = histcounts(AI,0:(parameters.AgeDeath));
    SummaryStatistics.AgeVaccIncidence(:,5, CurrentDay) = ...
        SummaryStatistics.AgeVaccIncidence(:,5, CurrentDay) + NumInfectedEachAgeGroup(:);
        
    SummaryStatistics.InfectorInfecteePairs = ...
        [SummaryStatistics.InfectorInfecteePairs; InfectorInfecteePairsThisTimeStep];
    SummaryStatistics.ClinicalPathwaysData = ...
        [SummaryStatistics.ClinicalPathwaysData; ClinicalPathwaysData];

end

% Update agents due to aging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AgentCharacteristics, parameters] = aging(parameters,AgentCharacteristics)

    AgentCharacteristics.Age = parameters.dt_years + ...
                    AgentCharacteristics.Age;

    % Death and Birth: remove agents that are older than
    % AgeDeath.  Replace with newborns. 
    
    D=[];       
    D = [D; find(AgentCharacteristics.Age > parameters.AgeDeath)];
    
    % Remove agents that die
    AgentCharacteristics.InfectionStatus(D,:) = [];
    AgentCharacteristics.VaccinationStatus(D,:) = [];
    AgentCharacteristics.Age(D) = [];
    AgentCharacteristics.HouseholdList(D,:) = [];
    AgentCharacteristics.CurrentHousehold(D) = [];
    AgentCharacteristics.ID(D) = [];
    AgentCharacteristics.Contacts(D,:) = [];
    AgentCharacteristics.TestSensitivity(D,:) = [];
    AgentCharacteristics.Symptoms(D,:) = [];
    AgentCharacteristics.Quarantine(D,:) = [];
    AgentCharacteristics.Lockdown(D,:) = [];
    AgentCharacteristics.ReactiveVaccination(D,:) = [];

    % Include agents that are born
    NumberBirths = length(D);
    if NumberBirths > 0
 
        % Determine their current house and household list
        newHouseholdList = zeros(NumberBirths,3);
        newCurrentHousehold = zeros(NumberBirths,1);
        newIDs = zeros(NumberBirths,1);

        for j = 1: NumberBirths
            
            % Random numbers for household selection
            rns = parameters.UniformRand(parameters.countUR:parameters.countUR+5-1,1);
            rncount = 1;
            
            % Update counter for random numbers
            parameters.countUR = parameters.countUR + 5;
            if parameters.countUR > 0.9*10^5
                parameters = updatecounter(parameters, "uniform");
            end
            
            % New unique agent IDs, current household and household list 
            % for newborns
            newIDs(j) = parameters.countID + j;           
            newhh = ceil(rns(rncount)*parameters.NumberHouses); % randomly select
            newCurrentHousehold(j) = newhh;
            hhl = ceil(rns(rncount+1:rncount+3)*parameters.NumberHouses); % randomly select
            newHouseholdList(j,:) = parameters.HHIDs(hhl)';
            
            % Update counters for number agents in each house
            parameters.nahh = histcounts(AgentCharacteristics.CurrentHousehold,1:parameters.NumberHouses+1);
        
        end
        
        % Store data for newborns (assume susceptible, not vaccinated, not 
        % in islation, not in quarantine, but in lockdown if community currenlty is in lockdown)
        AgentCharacteristics.InfectionStatus = [AgentCharacteristics.InfectionStatus; zeros(NumberBirths,4)];
        AgentCharacteristics.VaccinationStatus = [AgentCharacteristics.VaccinationStatus; zeros(NumberBirths,5)];
        AgentCharacteristics.Age = [AgentCharacteristics.Age; 0.001 * ones(NumberBirths,1)];
        AgentCharacteristics.HouseholdList = [AgentCharacteristics.HouseholdList; newHouseholdList];
        AgentCharacteristics.CurrentHousehold = [AgentCharacteristics.CurrentHousehold; newCurrentHousehold];
        AgentCharacteristics.ID = [AgentCharacteristics.ID; newIDs];
        AgentCharacteristics.Contacts{end+NumberBirths,1} = [];
        
        AgentCharacteristics.TestSensitivity = [AgentCharacteristics.TestSensitivity; zeros(NumberBirths,1)];
        AgentCharacteristics.Symptoms = [AgentCharacteristics.Symptoms; zeros(NumberBirths,3)];
        AgentCharacteristics.Quarantine = [AgentCharacteristics.Quarantine; zeros(NumberBirths,2)];
        
        Hesitancy = zeros(NumberBirths,1);
        Hesitancy(rand(NumberBirths,1)<parameters.Hesitancy) = 1;    
        AgentCharacteristics.ReactiveVaccination = [AgentCharacteristics.ReactiveVaccination; ...
            zeros(NumberBirths,1) -ones(NumberBirths,1) Hesitancy];
        
        % If lockdown is on, place births into lockdown
        if parameters.response == 1
            if parameters.lockdown == 1
                newLD = [5*ones(NumberBirths,1) 0*ones(NumberBirths,1)];
                AgentCharacteristics.Lockdown = [AgentCharacteristics.Lockdown; newLD];
            else
                AgentCharacteristics.Lockdown = [AgentCharacteristics.Lockdown; zeros(NumberBirths,2)];
            end
        else
            AgentCharacteristics.Lockdown = [AgentCharacteristics.Lockdown; zeros(NumberBirths,2)];
        end
        
        parameters.countID = max(AgentCharacteristics.ID);
        parameters.PopSize = length(AgentCharacteristics.ID);
        
    end

end

% Update population change due to migration between dwellings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters, AgentCharacteristics] = hhmobility(parameters,AgentCharacteristics, CurrentHousehold, CurrentSymptomStatus, CurrentQuarantine, CurrentLockdown)

    IAgentsAll = (1:1:length(AgentCharacteristics.ID))';
    NewCurrentHH = CurrentHousehold;

    RelevantAgents = IAgentsAll;
    
    % Determine which agents do not move from current household as they are
    % in self isolation, isolation, quarantine or lockdown
    SelfIsoAgents = (CurrentSymptomStatus(:,1)==1 | CurrentSymptomStatus(:,1)==2);
    IsoAgents = (CurrentSymptomStatus(:,1)>=3);
    QAgents = (CurrentQuarantine(:,1)>=2);
    LDAgents = (CurrentLockdown(:,1)>=2);
    
    AgentsNoMobility = (SelfIsoAgents | IsoAgents | QAgents | LDAgents);
    RelevantAgents(AgentsNoMobility) = [];
    
    NumberRelevantAgents = length(RelevantAgents);

    % Random numbers for mobility decisions
    rns = parameters.UniformRand(parameters.countUR:parameters.countUR+NumberRelevantAgents-1,1);

    % Update counters for pregenerated random numbers
    parameters.countUR = parameters.countUR + NumberRelevantAgents;
    if parameters.countUR > 0.9*10^5
        parameters = updatecounter(parameters, "uniform");
    end

    % For each agent, use random numbers rns to determine whether they
    % stay in core household, regular, on/off or a
    % random household.
    IAgentsCore = RelevantAgents(rns<parameters.HHMobilityPD(1));
    IagentsReg = RelevantAgents(~ismember(RelevantAgents,IAgentsCore) & rns<parameters.HHMobilityPD(2));
    IagentsSpor = RelevantAgents(rns>parameters.HHMobilityPD(3));
    IagentsOnoff = RelevantAgents(~ismember(RelevantAgents,[IAgentsCore; IagentsReg; IagentsSpor]));

    NewCurrentHH(IAgentsCore,1) = AgentCharacteristics.HouseholdList(IAgentsCore,1);
    NewCurrentHH(IagentsReg,1) = AgentCharacteristics.HouseholdList(IagentsReg,2);
    NewCurrentHH(IagentsOnoff,1) = AgentCharacteristics.HouseholdList(IagentsOnoff,3);
    SporadicHouses = ceil(parameters.UniformRand(parameters.countUR:parameters.countUR+length(IagentsSpor)-1,1)*parameters.NumberHouses);
    NewCurrentHH(IagentsSpor,1) = parameters.HHIDs(SporadicHouses);
    
    % Agents in Q/LD/I/SI stay in core household
    NewCurrentHH(AgentsNoMobility,1) = AgentCharacteristics.HouseholdList(AgentsNoMobility,1);

    % Update counter for random numbers
    parameters.countUR = parameters.countUR + length(IagentsSpor);
    if parameters.countUR > 0.9*10^5
        parameters = updatecounter(parameters, "uniform");
    end

    % Update current household
    AgentCharacteristics.CurrentHousehold = NewCurrentHH;

end


% Calculate summary statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SummaryStatistics = generate_summary_statistics(SummaryStatistics,timestep,parameters,AgentCharacteristics)
   
    for i=1:parameters.NumberHouses
        % Find current residents
        current_residents = (AgentCharacteristics.CurrentHousehold == i);
        % Of these, count up which ones are susceptible
        numsusc = sum(AgentCharacteristics.InfectionStatus(current_residents,1)==0);
        SummaryStatistics.NumberSusceptibleInDwellingsTime(i,timestep) = numsusc;   
    end
    
    % Tracking infection and vaccination status over time
    % For each infection state
    for i = 1:7
        % For each vaccination state >=1 dose
        for j = 1:4
            if j<=2
                k = 1; % Pfizer
                jj = j; % Dose #
            else
                k = 3; % AZ
                jj = j-2; % Dose #
            end
            % Get relevant infection state data for all agents
            IStemp = squeeze(AgentCharacteristics.InfectionStatus(:,1));
            % Get relevant vaccination state data for all agents
            VStemp = squeeze(AgentCharacteristics.VaccinationStatus(:,k));
            % Get age data for all agents
            Agestemp = AgentCharacteristics.Age;
            % Remove data of agents that don't satisfy relevant inf and 
            % vacc status
            if i<=6 % E, PI, I, PA, A, R
                Agestemp(IStemp~=i | VStemp~=jj)=[];
                %IStemp(IStemp~=i | VStemp~=jj)=[];
            else % S
                Agestemp(IStemp~=0 | VStemp~=jj)=[];
                %IStemp(IStemp~=0 | VStemp~=jj)=[];
            end
            % Count remaining agents
            SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(i,j,:,timestep) = histcounts(Agestemp,0:(parameters.AgeDeath));
        end
        % For no vaccine: repeat
        IStemp = squeeze(AgentCharacteristics.InfectionStatus(:,1));
        VStemp1 = squeeze(AgentCharacteristics.VaccinationStatus(:,1));
        VStemp2 = squeeze(AgentCharacteristics.VaccinationStatus(:,3));
        Agestemp = AgentCharacteristics.Age;
        if i<=6
            Agestemp(IStemp~=i | VStemp1~=0 | VStemp2~=0)=[];
            %IStemp(IStemp~=i | VStemp1~=0 | VStemp2~=0)=[];
        else
            Agestemp(IStemp~=0 | VStemp1~=0 | VStemp2~=0)=[];
            %IStemp(IStemp~=0 | VStemp1~=0 | VStemp2~=0)=[];
        end
        SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(i,5,:,timestep) = histcounts(Agestemp,0:(parameters.AgeDeath));
    end
    
    if mod((timestep - 1)*parameters.dt,1)==0
        SummaryStatistics.QuarantinePersonDays = [SummaryStatistics.QuarantinePersonDays;
                                                sum(AgentCharacteristics.Quarantine(:,1)>1)];
    end
    
    IsoAgents = sum(AgentCharacteristics.Symptoms(:,1)>=3);
    QAgents = sum(AgentCharacteristics.Quarantine(:,1)>=2);
    LDAgents = sum(AgentCharacteristics.Lockdown(:,1)>=2);
    
    SummaryStatistics.NumberInQTime(timestep,1) = (QAgents);
    SummaryStatistics.NumberInLTime(timestep,1) = (LDAgents);
    SummaryStatistics.NumberInITime(timestep,1) = (IsoAgents);    

end


% Pregenerate random numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parameters = prn(parameters)

    % Generate uniform random numbers for sampling contacts and mobility
    parameters.UniformRand = rand(1e5,1);
    parameters.countUR = 1;

    % Generate poisson random numbers for number of contacts
    parameters.ContactsNumberRand = zeros(parameters.NumberAgeClassesContacts,parameters.NumberAgeClassesContacts,1e5);
    for ii = 1: parameters.NumberAgeClassesContacts
        for jj = 1: parameters.NumberAgeClassesContacts
            parameters.ContactsNumberRand(ii,jj,:)=poissrnd(parameters.Ncontacts(ii,jj),[1e5,1]);
        end
    end     
    parameters.countCNR = 1;
    
    % Samppled duration of latency and incubation periods, calculated 
    % presymptomatic period, and sampled duration of symptomatic infection
    n = 1e5; 
    lp = lognrnd(parameters.MeanDurationLatency,parameters.SDLatency,n,1);
    ip = lognrnd(parameters.MeanDurationIncubation,parameters.SDIncubation,n,1);
    ln = lp(lp<ip); in = ip(lp<ip);
    s = lognrnd(parameters.MeanDuratioSymptoms,parameters.SDSymptoms,length(in),1);
    % Store infection duration sets: latency, pre-symp infectious, sympt
    % infectious
    parameters.DurationInfectionRand = [ln in-ln s];
    parameters.countDIR = 1;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameters = updatecounter(parameters, type)

    if type == "uniform"

        parameters.UniformRand = rand(1e5,1);
        parameters.countUR = 1;

    end
    
    if type == "durationinfection"
       
        n = 1e5; 
        % Sampled duration of latency and incubation periods, calculated 
        % presymptomatic period, and sampled duration of symptomatic infection
        lp = lognrnd(parameters.MeanDurationLatency,parameters.SDLatency,n,1);
        ip = lognrnd(parameters.MeanDurationIncubation,parameters.SDIncubation,n,1);
        ln = lp(lp<ip); in = ip(lp<ip);
        s = lognrnd(parameters.MeanDuratioSymptoms,parameters.SDSymptoms,length(in),1);
        % Store infection duration sets: latency, pre-symp infectious, sympt
        % infectious
        parameters.DurationInfectionRand = [ln in-ln s];
        parameters.countDIR = 1;
        
    end

    if type == "contacts"
        
        for ii = 1: parameters.NumberAgeClassesContacts
            for jj = 1: parameters.NumberAgeClassesContacts
                parameters.ContactsNumberRand(ii,jj,:)=poissrnd(parameters.Ncontacts(ii,jj),[1e5,1]);
            end
        end
        parameters.countCNR = 1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AgeGroup = calculate_age_groups(Ages,ACD)
    AgeGroup = zeros(size(Ages));
    for pp = 1: length(ACD)-1
        AgeGroup(Ages >= ACD(pp) & Ages < ACD(pp+1)) = pp;
    end
end

