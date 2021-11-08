function [AgentCharacteristics, SummaryStatistics, parameters] = ...
    manage_quarantine(parameters,AgentCharacteristics,SummaryStatistics,...
    CurrentSymptomStatus, CurrentTestSensitivity,CurrentContacts, CurrentInfection,...
    CurrentQuarantine, KeepIndexAllAgents, CurrentIDs, timestep, CurrentHousehold, CurrentAge, CurrentVaccinationStatus)

% Agents who enter quarantine are tested. If they receive a positive test, 
% their contacts are traced, and they enter isolation after a delay.  If 
% negative result, they remain in quarantine.  Before they are due to be 
% release from quarantine they are tested again.  If they receive a negative 
% test, they can leave quarantine, otherwise they enter isolation.

KeepIndexAllAgents = KeepIndexAllAgents(:);

% Determine agent who are still waiting entry into quarantine
IndexAgentsWaiting = KeepIndexAllAgents(CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) > 0);

%Determine agents who enter quarantine, and their entry test result (this 
% requires specifying time since (or before) symptom onset at time of test
Susc = KeepIndexAllAgents(CurrentInfection(:,1)==0 & CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) <= 0);
STSSO = zeros(size(Susc));
Exposed = KeepIndexAllAgents(CurrentInfection(:,1)==1 & CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) <= 0);
ETSSO = -(CurrentInfection(Exposed,2) + CurrentInfection(Exposed,3));
PreI = KeepIndexAllAgents(CurrentInfection(:,1)==2 & CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) <= 0);
PITSSO = -CurrentInfection(PreI,3);
I = KeepIndexAllAgents(CurrentInfection(:,1)==3 & CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) <= 0);
ITSSO = CurrentSymptomStatus(I,3);
PreA = KeepIndexAllAgents(CurrentInfection(:,1)==4 & CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) <= 0);
PATSSO = -CurrentInfection(PreA,3);
A = KeepIndexAllAgents(CurrentInfection(:,1)==5 & CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) <= 0);
ATSSO = CurrentSymptomStatus(A,3);
R = KeepIndexAllAgents(CurrentInfection(:,1)==6 & CurrentQuarantine(:,1)==1 & CurrentQuarantine(:,2) <= 0);
RTSSO = CurrentSymptomStatus(R,3);

IndexAgentsEnter = [Susc(:); Exposed(:); PreI(:); I(:); PreA(:); A(:); R(:)];
TSSO = [STSSO(:); ETSSO(:); PITSSO(:); ITSSO(:); PATSSO(:); ATSSO(:); RTSSO(:)]; 

[TestSensitivityAwaitingTest, parameters] = test_sensitivity(TSSO,parameters);

TestSensitivityAwaitingTest = TestSensitivityAwaitingTest .* CurrentTestSensitivity(IndexAgentsEnter);

RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexAgentsEnter)-1);

% Negative result
IndexAgentsFullQuarantine = IndexAgentsEnter(RandomNumbers > TestSensitivityAwaitingTest);
% True positive
IndexAgentsEnterIso = IndexAgentsEnter(RandomNumbers <= TestSensitivityAwaitingTest);

% Update counter for pregenerated random numbers:
parameters.countUR = parameters.countUR + length(IndexAgentsEnter);
if parameters.countUR > 0.9*10^5
    parameters.UniformRand = rand(1e5,1);
    parameters.countUR = 1;
end

% Update AC & SS
AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) = AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) - parameters.dt;
AgentCharacteristics.Quarantine(IndexAgentsFullQuarantine,1) = 4;
AgentCharacteristics.Quarantine(IndexAgentsFullQuarantine,2) = parameters.DelayOnsetSymptTestResult;
AgentCharacteristics.Quarantine(IndexAgentsEnterIso,1) = 2;
AgentCharacteristics.Quarantine(IndexAgentsEnterIso,2) = parameters.DelayOnsetSymptTestResult;

SummaryStatistics.TestDate = [SummaryStatistics.TestDate;
    CurrentIDs(IndexAgentsEnterIso) ones(size(IndexAgentsEnterIso))*(timestep - 1)*parameters.dt];

% Determine agents who are still awaiting entry test result
% Those heading to isolation:
IndexAgentsWaiting1 = KeepIndexAllAgents(CurrentQuarantine(:,1)==2 & CurrentQuarantine(:,2) > 0);
% Those who will be confirmed to stay in quarantine
IndexAgentsWaiting2 = KeepIndexAllAgents(CurrentQuarantine(:,1)==4 & CurrentQuarantine(:,2) > 0);

% Determine agents who receive entry test result 
IndexAgentsGoIso = KeepIndexAllAgents(CurrentQuarantine(:,1)==2 & CurrentQuarantine(:,2) <= 0);
IndexAgentsConfirmed = KeepIndexAllAgents(CurrentQuarantine(:,1)==4 & CurrentQuarantine(:,2) <= 0);

% Find and track contacts of +ve agents (immediately on knowing test result)
if ~isempty(IndexAgentsGoIso)
    [AgentCharacteristics, parameters] = ...
        trace_contacts(parameters, AgentCharacteristics, CurrentContacts, CurrentSymptomStatus, ...
        CurrentInfection, IndexAgentsGoIso, KeepIndexAllAgents, CurrentIDs, timestep-1);
end

% Update AC
AgentCharacteristics.Quarantine(IndexAgentsWaiting1,2) = AgentCharacteristics.Quarantine(IndexAgentsWaiting1,2) - parameters.dt;
AgentCharacteristics.Quarantine(IndexAgentsWaiting2,2) = AgentCharacteristics.Quarantine(IndexAgentsWaiting2,2) - parameters.dt;
AgentCharacteristics.Quarantine(IndexAgentsGoIso,1) = 3;
AgentCharacteristics.Quarantine(IndexAgentsGoIso,2) = parameters.DelayTestResultIsolation;
AgentCharacteristics.Quarantine(IndexAgentsConfirmed,1) = 5;
AgentCharacteristics.Quarantine(IndexAgentsConfirmed,2) = parameters.QuarantineClearanceTestDay - parameters.DelayOnsetSymptTestResult;

% Determine agents who still waiting to go to isolation
IndexAgentsWaiting = KeepIndexAllAgents(CurrentQuarantine(:,1)==3 & CurrentQuarantine(:,2) > 0);
% Determine agents who go to isolation
IndexAgentsGoIso = KeepIndexAllAgents(CurrentQuarantine(:,1)==3 & CurrentQuarantine(:,2) <= 0);

% Update AC
AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) = AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) - parameters.dt;
%For agents going to isolation, remove from quarantine pathway, add to isolation pathway
AgentCharacteristics.Quarantine(IndexAgentsGoIso,1) = 0;
AgentCharacteristics.Quarantine(IndexAgentsGoIso,2) = 0;
AgentCharacteristics.Symptoms(IndexAgentsGoIso,1) = 3;
AgentCharacteristics.Symptoms(IndexAgentsGoIso,2) = parameters.IsolationClearanceTestDay;

newdata = [CurrentIDs(IndexAgentsGoIso) CurrentHousehold(IndexAgentsGoIso) CurrentAge(IndexAgentsGoIso)...
                                (timestep-1)*parameters.dt*ones(size(IndexAgentsGoIso)) 2*ones(size(IndexAgentsGoIso))...
                                CurrentVaccinationStatus(IndexAgentsGoIso,:)];

SummaryStatistics.CaseData = [SummaryStatistics.CaseData; newdata];

% Determine agents who are yet to reach clearance test day
IndexAgentsWaiting = KeepIndexAllAgents(CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) > 0);

% Determine agents who reach clearance test day, and results of their test
Susc = KeepIndexAllAgents(CurrentInfection(:,1)==0 & CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) <= 0);
STSSO = zeros(size(Susc));
Exposed = KeepIndexAllAgents(CurrentInfection(:,1)==1 & CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) <= 0);
ETSSO = -(CurrentInfection(Exposed,2) + CurrentInfection(Exposed,3));
PreI = KeepIndexAllAgents(CurrentInfection(:,1)==2 & CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) <= 0);
PITSSO = -CurrentInfection(PreI,3);
I = KeepIndexAllAgents(CurrentInfection(:,1)==3 & CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) <= 0);
ITSSO = CurrentSymptomStatus(I,3);
PreA = KeepIndexAllAgents(CurrentInfection(:,1)==4 & CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) <= 0);
PATSSO = -CurrentInfection(PreA,3);
A = KeepIndexAllAgents(CurrentInfection(:,1)==5 & CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) <= 0);
ATSSO = CurrentSymptomStatus(A,3);
R = KeepIndexAllAgents(CurrentInfection(:,1)==6 & CurrentQuarantine(:,1)==5 & CurrentQuarantine(:,2) <= 0);
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
AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) = AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) - parameters.dt;
AgentCharacteristics.Quarantine(IndexAgentsToLeave,1) = 6;
AgentCharacteristics.Quarantine(IndexAgentsToLeave,2) = parameters.MinDurationQua - parameters.QuarantineClearanceTestDay;
AgentCharacteristics.Quarantine(IndexAgentsToEnterIso,1) = 2;
AgentCharacteristics.Quarantine(IndexAgentsToEnterIso,2) = parameters.DelayOnsetSymptTestResult;

SummaryStatistics.TestDate = [SummaryStatistics.TestDate;
    CurrentIDs(IndexAgentsToEnterIso) ones(size(IndexAgentsToEnterIso))*(timestep - 1)*parameters.dt];

% Find agents who are waiting on release, and who leave quarantine at 14 days
IndexAgentsWaiting = KeepIndexAllAgents(CurrentQuarantine(:,1)==6 & CurrentQuarantine(:,2) > 0);
IndexAgentsLeaving = KeepIndexAllAgents(CurrentQuarantine(:,1)==6 & CurrentQuarantine(:,2) <= 0);

% Update AC
AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) = AgentCharacteristics.Quarantine(IndexAgentsWaiting,2) - parameters.dt;
AgentCharacteristics.Quarantine(IndexAgentsLeaving,1) = 0;
AgentCharacteristics.Quarantine(IndexAgentsLeaving,2) = 0;

