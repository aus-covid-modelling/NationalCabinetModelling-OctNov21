function [AgentCharacteristics, SummaryStatistics, parameters] = ...
    lockdown(parameters,AgentCharacteristics,SummaryStatistics,...
    CurrentSymptomStatus, CurrentTestSensitivity,CurrentContacts, CurrentInfection,...
    CurrentLockdown, KeepIndexAllAgents, CurrentIDs, timestep, CurrentHousehold, CurrentAge, CurrentVaccinationStatus)

% Agents who enter lockdown are tested. If they receive a positive test, 
% their contacts are traced, and they enter isolation after a delay.  If 
% negative result, they remain in lockdown.  Before they are due to be 
% release from lockdown they are tested again.  If they receive a negative 
% test, they can leave lockdown, otherwise they enter isolation.

KeepIndexAllAgents = KeepIndexAllAgents(:);

% Determine agents who are still waiting entry into lockdown
IndexAgentsWaiting = KeepIndexAllAgents(CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) > 0);

% Determine agents who enter lockdown, and their entry test result (this 
% requires specifying time since (or before) symptom onset at time of test
Susc = KeepIndexAllAgents(CurrentInfection(:,1)==0 & CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) <= 0);
STSSO = zeros(size(Susc));
Exposed = KeepIndexAllAgents(CurrentInfection(:,1)==1 & CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) <= 0);
ETSSO = -(CurrentInfection(Exposed,2) + CurrentInfection(Exposed,3));
PreI = KeepIndexAllAgents(CurrentInfection(:,1)==2 & CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) <= 0);
PITSSO = -CurrentInfection(PreI,3);
I = KeepIndexAllAgents(CurrentInfection(:,1)==3 & CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) <= 0);
ITSSO = CurrentSymptomStatus(I,3);
PreA = KeepIndexAllAgents(CurrentInfection(:,1)==4 & CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) <= 0);
PATSSO = -CurrentInfection(PreA,3);
A = KeepIndexAllAgents(CurrentInfection(:,1)==5 & CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) <= 0);
ATSSO = CurrentSymptomStatus(A,3);
R = KeepIndexAllAgents(CurrentInfection(:,1)==6 & CurrentLockdown(:,1)==1 & CurrentLockdown(:,2) <= 0);
RTSSO = CurrentSymptomStatus(R,3);

IndexAgentsEnter = [Susc(:); Exposed(:); PreI(:); I(:); PreA(:); A(:); R(:)];
TSSO = [STSSO(:); ETSSO(:); PITSSO(:); ITSSO(:); PATSSO(:); ATSSO(:); RTSSO(:)]; 

% Determine test sensitivity based on time since (or before) symptom onset
[TestSensitivityAwaitingTest, parameters] = test_sensitivity(TSSO,parameters);

% Modify this for susceptible and recovered (where test senstity is zero) 
TestSensitivityAwaitingTest = TestSensitivityAwaitingTest .* CurrentTestSensitivity(IndexAgentsEnter);

% Now determine test result
RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexAgentsEnter)-1);
% Negative result
IndexAgentsFullLockdown = IndexAgentsEnter(RandomNumbers > TestSensitivityAwaitingTest);
% True positive
IndexAgentsEnterIso = IndexAgentsEnter(RandomNumbers <= TestSensitivityAwaitingTest);

% Update counter for pregenerated random numbers:
parameters.countUR = parameters.countUR + length(IndexAgentsEnter);
if parameters.countUR > 0.9*10^5
    parameters.UniformRand = rand(1e5,1);
    parameters.countUR = 1;
end

% Update AC & SS
AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) = AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) - parameters.dt;
AgentCharacteristics.Lockdown(IndexAgentsFullLockdown,1) = 4;
AgentCharacteristics.Lockdown(IndexAgentsFullLockdown,2) = parameters.DelayLockdownTestResultsCommunity;
AgentCharacteristics.Lockdown(IndexAgentsEnterIso,1) = 2;
AgentCharacteristics.Lockdown(IndexAgentsEnterIso,2) = parameters.DelayLockdownTestResultsCommunity;

SummaryStatistics.TestDate = [SummaryStatistics.TestDate;
    CurrentIDs(IndexAgentsEnterIso) ones(size(IndexAgentsEnterIso))*(timestep - 1)*parameters.dt];

% Determine agents who are still awaiting entry test result
% Those heading to isolation:
IndexAgentsWaiting1 = KeepIndexAllAgents(CurrentLockdown(:,1)==2 & CurrentLockdown(:,2) > 0);
% Those who will be confirmed to stay in lockdown
IndexAgentsWaiting2 = KeepIndexAllAgents(CurrentLockdown(:,1)==4 & CurrentLockdown(:,2) > 0);

% Determine agents who receive entry test result 
IndexAgentsGoIso = KeepIndexAllAgents(CurrentLockdown(:,1)==2 & CurrentLockdown(:,2) <= 0);
IndexAgentsConfirmed = KeepIndexAllAgents(CurrentLockdown(:,1)==4 & CurrentLockdown(:,2) <= 0);

% Find and track contacts of +ve agents (immediately on knowing test result)
if ~isempty(IndexAgentsGoIso)
    [AgentCharacteristics, parameters] = ...
        trace_contacts(parameters, AgentCharacteristics, CurrentContacts, CurrentSymptomStatus, ...
        CurrentInfection, IndexAgentsGoIso, KeepIndexAllAgents, CurrentIDs, timestep-1);
end

% Update AC
AgentCharacteristics.Lockdown(IndexAgentsWaiting1,2) = AgentCharacteristics.Lockdown(IndexAgentsWaiting1,2) - parameters.dt;
AgentCharacteristics.Lockdown(IndexAgentsWaiting2,2) = AgentCharacteristics.Lockdown(IndexAgentsWaiting2,2) - parameters.dt;
AgentCharacteristics.Lockdown(IndexAgentsGoIso,1) = 3;
AgentCharacteristics.Lockdown(IndexAgentsGoIso,2) = parameters.DelayTestResultIsolation;
AgentCharacteristics.Lockdown(IndexAgentsConfirmed,1) = 5;
AgentCharacteristics.Lockdown(IndexAgentsConfirmed,2) = parameters.LockdownClearanceTestDay - parameters.DelayLockdownTestResultsCommunity;

% Determine agents who still waiting to go to isolation
IndexAgentsWaiting = KeepIndexAllAgents(CurrentLockdown(:,1)==3 & CurrentLockdown(:,2) > 0);
% Determine agents who go to isolation
IndexAgentsGoIso = KeepIndexAllAgents(CurrentLockdown(:,1)==3 & CurrentLockdown(:,2) <= 0);

% Update AC
AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) = AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) - parameters.dt;
% For agents going to isolation, remove from quarantine pathway, add to isolation
% pathway but keep lockdown details
AgentCharacteristics.Quarantine(IndexAgentsGoIso,1) = 0;
AgentCharacteristics.Quarantine(IndexAgentsGoIso,2) = 0;
AgentCharacteristics.Symptoms(IndexAgentsGoIso,1) = 3;
AgentCharacteristics.Symptoms(IndexAgentsGoIso,2) = parameters.IsolationClearanceTestDay;
AgentCharacteristics.Lockdown(IndexAgentsGoIso,1) = 5;
AgentCharacteristics.Lockdown(IndexAgentsGoIso,2) = parameters.LockdownClearanceTestDay - ...
    parameters.DelayLockdownTestResultsCommunity - parameters.DelayTestResultIsolation;

% Update case data
SummaryStatistics.CaseData = [SummaryStatistics.CaseData;
    CurrentIDs(IndexAgentsGoIso) CurrentHousehold(IndexAgentsGoIso) CurrentAge(IndexAgentsGoIso)...
                                (timestep-1)*parameters.dt*ones(size(IndexAgentsGoIso)) 3*ones(size(IndexAgentsGoIso))...
                                CurrentVaccinationStatus(IndexAgentsGoIso,:)];

% Determine agents who are yet to reach clearance test day
IndexAgentsWaiting = KeepIndexAllAgents(CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) > 0);

% Determine agents who reach clearance test day, and results of their test
Susc = KeepIndexAllAgents(CurrentInfection(:,1)==0 & CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) <= 0);
STSSO = zeros(size(Susc));
Exposed = KeepIndexAllAgents(CurrentInfection(:,1)==1 & CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) <= 0);
ETSSO = -(CurrentInfection(Exposed,2) + CurrentInfection(Exposed,3));
PreI = KeepIndexAllAgents(CurrentInfection(:,1)==2 & CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) <= 0);
PITSSO = -CurrentInfection(PreI,3);
I = KeepIndexAllAgents(CurrentInfection(:,1)==3 & CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) <= 0);
ITSSO = CurrentSymptomStatus(I,3);
PreA = KeepIndexAllAgents(CurrentInfection(:,1)==4 & CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) <= 0);
PATSSO = -CurrentInfection(PreA,3);
A = KeepIndexAllAgents(CurrentInfection(:,1)==5 & CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) <= 0);
ATSSO = CurrentSymptomStatus(A,3);
R = KeepIndexAllAgents(CurrentInfection(:,1)==6 & CurrentLockdown(:,1)==5 & CurrentLockdown(:,2) <= 0);
RTSSO = CurrentSymptomStatus(R,3);

IndexAgentsTested = [Susc(:); Exposed(:); PreI(:); I(:); PreA(:); A(:); R(:)];
TSSO = [STSSO(:); ETSSO(:); PITSSO(:); ITSSO(:); PATSSO(:); ATSSO(:); RTSSO(:)]; 

[TestSensitivityAwaitingTest, parameters] = test_sensitivity(TSSO,parameters);
TestSensitivityAwaitingTest = TestSensitivityAwaitingTest .* CurrentTestSensitivity(IndexAgentsTested);

RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexAgentsTested)-1);

% Negative result
IndexAgentsToLeave = IndexAgentsTested(RandomNumbers > TestSensitivityAwaitingTest);
% True positive
IndexAgentsToEnterIso = IndexAgentsTested(RandomNumbers <= TestSensitivityAwaitingTest);

% Update counter for pregenerated random numbers:
parameters.countUR = parameters.countUR + length(IndexAgentsEnter);
if parameters.countUR > 0.9*10^5
    parameters.UniformRand = rand(1e5,1);
    parameters.countUR = 1;
end

% Update AC
AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) = AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) - parameters.dt;
AgentCharacteristics.Lockdown(IndexAgentsToLeave,1) = 6;
AgentCharacteristics.Lockdown(IndexAgentsToLeave,2) = parameters.MinDurationLoc - parameters.LockdownClearanceTestDay;
AgentCharacteristics.Lockdown(IndexAgentsToEnterIso,1) = 2;
AgentCharacteristics.Lockdown(IndexAgentsToEnterIso,2) = parameters.DelayLockdownTestResultsCommunity;

SummaryStatistics.TestDate = [SummaryStatistics.TestDate;
    CurrentIDs(IndexAgentsToEnterIso) ones(size(IndexAgentsToEnterIso))*(timestep - 1)*parameters.dt];

% Don't remove agents from lockdown here, this is controlled in simulator.m
IndexAgentsWaiting = KeepIndexAllAgents(CurrentLockdown(:,1)==6 & CurrentLockdown(:,2) > 0);
% Update AC
AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) = AgentCharacteristics.Lockdown(IndexAgentsWaiting,2) - parameters.dt;




