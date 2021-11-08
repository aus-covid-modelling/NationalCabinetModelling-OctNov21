% This script will run the model once according to parameters set in:
% parameters_base
% parameters_infection (with exception of intial vacc status which are overwritten here)
% parameters_response

parameters=[]; 

% Choose population
% 0: Size 1000, mean HH size 7.7 
% 1: Size 220, mean HH size 6.1 
% 2: Size 580, mean HH size 4.8 
% 3: Size 1018, mean HH size 3.5 
parameters.exemplar = 0; 

% Load population parameters
parameters = parameters_base(parameters);

% Load infection-related parameters
parameters = parameters_infection(parameters);

% Initial vaccintation proportions in age groups:
% Exemplars 0:2, 5 age groups: [0 12 16 40 60 parameters.AgeDeath]
% Exemplar 3, 10 age groups: [0 12 16 20:10:80 parameters.AgeDeath]
parameters.VP10 = [0 0.25 0.25 0.25 0]; % these should be length 10 for exemplar 3, 5 otherwise
parameters.VP20 = [0 0.25 0.25 0.25 0]; 
parameters.VAZ10 = [0 0 0 0 0.1]; 
parameters.VAZ20 = [0 0 0 0 0.7]; 

% Load response-related parameters
parameters = parameters_response(parameters);

% Initialise age structure and household structure
[AgentCharacteristics, parameters] = initialise_demographics(parameters);

% Initialise infection and immunity status
AgentCharacteristics = initialise_infection_status(AgentCharacteristics,parameters);

% Run simulation
[AgentCharacteristics, SummaryStatistics, parameters] = ...
                simulator(AgentCharacteristics, parameters);

% Integrate over age and/or vaccination status to get cumulative totals
% over time
TotalInfectiousTime = sum(squeeze(sum(squeeze(sum(SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(1:5,:,:,:),2)),1)),1);
TotalVaccinatedTime = sum(squeeze(sum(squeeze(sum(SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(:,1:4,:,:),2)),1)),1);
TotalVaccinated1DoseTime = sum(squeeze(sum(squeeze(sum(SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(:,[1 3],:,:),2)),1)),1);
TotalVaccinated2DoseTime = sum(squeeze(sum(squeeze(sum(SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(:,[2 4],:,:),2)),1)),1);
TotalRecoveredNonVacced = sum(squeeze(sum(squeeze(SummaryStatistics.NumberInfectionVaccinationAgeStatusTime(6,5,:,:)),1)),1);
TotalNotRecoveredNonVacced = sum(squeeze(sum(squeeze(SummaryStatistics.NumberInfectionVaccinationAgeStatusTime([1:5 7],5,:,:)),1)),1);

% Plot outputs for one simulation
figure
plot(TotalInfectiousTime/parameters.PopSize,'LineWidth',2)
hold on
plot(SummaryStatistics.NumberInQTime/parameters.PopSize,'LineWidth',2)
hold on
plot(SummaryStatistics.NumberInITime/parameters.PopSize,'LineWidth',2)
hold on
plot(SummaryStatistics.NumberInLTime/parameters.PopSize,'LineWidth',2)
hold on
plot(TotalVaccinatedTime/parameters.PopSize,'LineWidth',2)
hold on
plot(TotalVaccinated1DoseTime/parameters.PopSize,'LineWidth',2)
hold on
plot(TotalVaccinated2DoseTime/parameters.PopSize,'LineWidth',2)
hold on

legend('Infected','Quarantined','Isolated','In lockdown','Vaccinated', '1 Dose','2 Dose')
xlim([0 100])
ylim([0 1.1])
