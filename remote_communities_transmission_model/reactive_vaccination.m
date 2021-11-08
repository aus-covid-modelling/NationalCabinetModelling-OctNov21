function [AgentCharacteristics, SummaryStatistics, parameters] = ...
        reactive_vaccination(parameters, AgentCharacteristics, ...
                SummaryStatistics,CurrentInfection, CurrentVaccination, ...
                KeepIndexAllAgents,CurrentReactiveVaccination)

% This script simulates reactive vaccination (agents receiving first dose)            
            
% Find non-vaccinated agents who are currently suceptible, not vaccine
% hesitant, and old enough to be vaccinated
NonVacced = KeepIndexAllAgents(CurrentVaccination(:,1)==0 & ...
    CurrentVaccination(:,3)==0 & CurrentInfection(:,1)==0 & ...
    CurrentReactiveVaccination(:,3)==0 & ...
    AgentCharacteristics.Age>parameters.MinAgeVaccination);

% And their ages
AgesNonVacced = AgentCharacteristics.Age(NonVacced);

% If older first, sort by age
if parameters.OlderFirst == 1
    [~,ind] = sort(AgesNonVacced,'descend');
else
    % Otherwise randomly order
    ind = randperm(size(NonVacced,1));   
end
NonVacced = NonVacced(ind);
AgesNonVacced = AgesNonVacced(ind);

if ~isempty(NonVacced)
    % Select those to be vaccinated today
    if length(NonVacced) >= parameters.DailyMaxVaccination
        NonVacced = NonVacced(1:parameters.DailyMaxVaccination);
        AgesNonVacced = AgesNonVacced(1:parameters.DailyMaxVaccination);
    end
    
    % Give first dose
    % Pfizer
    AgentCharacteristics.VaccinationStatus(NonVacced(AgesNonVacced<60),1) = 1;
    AgentCharacteristics.VaccinationStatus(NonVacced(AgesNonVacced<60),2) = 0;
    % AZ
    AgentCharacteristics.VaccinationStatus(NonVacced(AgesNonVacced>=60),3) = 1;
    AgentCharacteristics.VaccinationStatus(NonVacced(AgesNonVacced>=60),4) = 0;

    % Update second dose timers
    AgentCharacteristics.ReactiveVaccination(NonVacced(AgesNonVacced<60),1) = parameters.DoseIntervalP;
    AgentCharacteristics.ReactiveVaccination(NonVacced(AgesNonVacced>=60),1) = parameters.DoseIntervalAZ;

end
