function [AgentCharacteristics, SummaryStatistics, parameters] = ...
    manage_agents_with_symptoms(parameters,AgentCharacteristics,SummaryStatistics,...
    CurrentSymptomStatus,CurrentTestSensitivity,CurrentContacts,CurrentInfection, ...
    KeepIndexAllAgents, CurrentIDs, timestep, CurrentHousehold, CurrentAge,CurrentVacc)

% Agents who develop symptoms self isolate until they obtain their test
% result.  If they receive a positive test, their contacts are traced, and 
% they enter isolation after a delay.  If negative result, they leave self 
% isolation.  Cases are tested before they are due to be release from
% isolation.  If they receive a negative test, they can leave isolation,
% otherwise they enter another round of isolation.

KeepIndexAllAgents = KeepIndexAllAgents(:);

% Determine agents who will stay in self isolation 
IndexAgentsStaySelfIso = KeepIndexAllAgents(CurrentSymptomStatus(:,1)==1 & CurrentSymptomStatus(:,2) > 0);

% Determine agents who will leave self isolation and their entry test result (this 
% requires specifying time since symptom onset at time of test
IndexAgentsLeaveSelfIso = KeepIndexAllAgents(CurrentSymptomStatus(:,1)==1 & CurrentSymptomStatus(:,2) <= 0);
TSSO = zeros(size(IndexAgentsLeaveSelfIso)); % time since symptom onset (assume got tested when symptoms developed)
[TestSensitivityLeaveIso, parameters] = test_sensitivity(TSSO,parameters);

TestSensitivityLeaveIso = TestSensitivityLeaveIso.* CurrentTestSensitivity(IndexAgentsLeaveSelfIso); % dont need this, except if want to turn isolation off easily

% Determine where they go by test result
RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexAgentsLeaveSelfIso)-1);
% Negative result
IndexAgentsLeaveSelfIsoNegResult = IndexAgentsLeaveSelfIso(RandomNumbers > TestSensitivityLeaveIso);
% True positive
IndexAgentsLeaveSelfIsoIso = IndexAgentsLeaveSelfIso(RandomNumbers <= TestSensitivityLeaveIso);

% Update counter for pregenerated random numbers:
parameters.countUR = parameters.countUR + length(IndexAgentsLeaveSelfIso);
if parameters.countUR > 0.9*10^5
    parameters.UniformRand = rand(1e5,1);
    parameters.countUR = 1;
end

% Update AC
AgentCharacteristics.Symptoms(IndexAgentsStaySelfIso,2) = AgentCharacteristics.Symptoms(IndexAgentsStaySelfIso,2) - parameters.dt;
AgentCharacteristics.Symptoms(IndexAgentsLeaveSelfIsoNegResult,1) = 0;
AgentCharacteristics.Symptoms(IndexAgentsLeaveSelfIsoNegResult,2) = 0;
AgentCharacteristics.Symptoms(IndexAgentsLeaveSelfIsoIso,1) = 2;
AgentCharacteristics.Symptoms(IndexAgentsLeaveSelfIsoIso,2) = parameters.DelayTestResultIsolation;

SummaryStatistics.TestDate = [SummaryStatistics.TestDate;
    CurrentIDs(IndexAgentsLeaveSelfIsoIso) ones(size(IndexAgentsLeaveSelfIsoIso))*(timestep - 1)*parameters.dt];

% Find and track contacts of +ve agents (immediately on knowing test result)
if ~isempty(IndexAgentsLeaveSelfIsoIso)
    [AgentCharacteristics, parameters] = ...
        trace_contacts(parameters, AgentCharacteristics, CurrentContacts, CurrentSymptomStatus, ...
        CurrentInfection, IndexAgentsLeaveSelfIsoIso, KeepIndexAllAgents, CurrentIDs, timestep-1);
end

% Determine agents who get moved from self isolation to isolation and update status
IndexAgentsLeaveSelfIso = KeepIndexAllAgents(CurrentSymptomStatus(:,1)==2 & CurrentSymptomStatus(:,2) <= 0);
IndexAgentsDontLeaveSelfIso = KeepIndexAllAgents(CurrentSymptomStatus(:,1)==2 & CurrentSymptomStatus(:,2) > 0);

% Update AC and SS
AgentCharacteristics.Symptoms(IndexAgentsDontLeaveSelfIso,2) = AgentCharacteristics.Symptoms(IndexAgentsDontLeaveSelfIso,2) - parameters.dt;
AgentCharacteristics.Symptoms(IndexAgentsLeaveSelfIso,1) = 3;
AgentCharacteristics.Symptoms(IndexAgentsLeaveSelfIso,2) = parameters.IsolationClearanceTestDay;

% For agents going to isolation, drop from quarantine pathway
% Before, find agents that are in Q for storing case data:
QS = AgentCharacteristics.Quarantine(IndexAgentsLeaveSelfIso,1);
IndexAgentsLeaveSelfIsoQ = IndexAgentsLeaveSelfIso(QS>=2);
IndexAgentsLeaveSelfIsoNotQ = IndexAgentsLeaveSelfIso(QS<2);

AgentCharacteristics.Quarantine(IndexAgentsLeaveSelfIso,1) = 0;
AgentCharacteristics.Quarantine(IndexAgentsLeaveSelfIso,2) = 0;

SummaryStatistics.CaseData = [SummaryStatistics.CaseData;
    CurrentIDs(IndexAgentsLeaveSelfIsoNotQ) CurrentHousehold(IndexAgentsLeaveSelfIsoNotQ) CurrentAge(IndexAgentsLeaveSelfIsoNotQ)...
                                (timestep-1)*parameters.dt*ones(size(IndexAgentsLeaveSelfIsoNotQ)) 1*ones(size(IndexAgentsLeaveSelfIsoNotQ))...
                                CurrentVacc(IndexAgentsLeaveSelfIsoNotQ,:)];
                            
SummaryStatistics.CaseData = [SummaryStatistics.CaseData;
    CurrentIDs(IndexAgentsLeaveSelfIsoQ) CurrentHousehold(IndexAgentsLeaveSelfIsoQ) CurrentAge(IndexAgentsLeaveSelfIsoQ)...
                                (timestep-1)*parameters.dt*ones(size(IndexAgentsLeaveSelfIsoQ)) 2*ones(size(IndexAgentsLeaveSelfIsoQ))...
                                CurrentVacc(IndexAgentsLeaveSelfIsoQ,:)];

% Determine agents who are yet to reach clearance test day
IndexAgentsNotAwaitingTest = (CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) > 0);

% Determine agents who reach clearance test day, who will transition to awaiting 
% their test reult, and the result of their test (this 
% requires specifying time since (or before) symptom onset at time of test
% TSSO calculations
Susc = KeepIndexAllAgents(CurrentInfection(:,1)==0 & CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) <= 0);
STSSO = zeros(size(Susc));
Exposed = KeepIndexAllAgents(CurrentInfection(:,1)==1 & CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) <= 0);
ETSSO = -(CurrentInfection(Exposed,2) + CurrentInfection(Exposed,3));
PreI = KeepIndexAllAgents(CurrentInfection(:,1)==2 & CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) <= 0);
PITSSO = -CurrentInfection(PreI,3);
I = KeepIndexAllAgents(CurrentInfection(:,1)==3 & CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) <= 0);
ITSSO = CurrentSymptomStatus(I,3);
PreA = KeepIndexAllAgents(CurrentInfection(:,1)==4 & CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) <= 0);
PATSSO = -CurrentInfection(PreA,3);
A = KeepIndexAllAgents(CurrentInfection(:,1)==5 & CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) <= 0);
ATSSO = CurrentSymptomStatus(A,3);
R = KeepIndexAllAgents(CurrentInfection(:,1)==6 & CurrentSymptomStatus(:,1)==3 & CurrentSymptomStatus(:,2) <= 0);
RTSSO = CurrentSymptomStatus(R,3);

IndexAgentsAwaitingTest = [Susc(:); Exposed(:); PreI(:); I(:); PreA(:); A(:); R(:)];
TSSO = [STSSO(:); ETSSO(:); PITSSO(:); ITSSO(:); PATSSO(:); ATSSO(:); RTSSO(:)]; 
[TestSensitivityAwaitingTest, parameters] = test_sensitivity(TSSO,parameters);
TestSensitivityAwaitingTest = TestSensitivityAwaitingTest .* CurrentTestSensitivity(IndexAgentsAwaitingTest);

RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(IndexAgentsAwaitingTest)-1);
% Negative result
IndexAgentsLeaveIso = IndexAgentsAwaitingTest(RandomNumbers > TestSensitivityAwaitingTest);
% True positive
IndexAgentsRestart = IndexAgentsAwaitingTest(RandomNumbers <= TestSensitivityAwaitingTest);

% Update counter for pregenerated random numbers:
parameters.countUR = parameters.countUR + length(IndexAgentsAwaitingTest);
if parameters.countUR > 0.9*10^5
    parameters.UniformRand = rand(1e5,1);
    parameters.countUR = 1;
end

% Update AC
AgentCharacteristics.Symptoms(IndexAgentsNotAwaitingTest,2) = AgentCharacteristics.Symptoms(IndexAgentsNotAwaitingTest,2) - parameters.dt;
AgentCharacteristics.Symptoms(IndexAgentsLeaveIso,1) = 4;
AgentCharacteristics.Symptoms(IndexAgentsLeaveIso,2) = parameters.MinDurationIso - parameters.IsolationClearanceTestDay;
AgentCharacteristics.Symptoms(IndexAgentsRestart,1) = 5;
AgentCharacteristics.Symptoms(IndexAgentsRestart,2) = parameters.DelayOnsetSymptTestResult;

% Determine agents who leave isolation 
IndexAgentsLeaveIso = (CurrentSymptomStatus(:,1)==4 & CurrentSymptomStatus(:,2) <= 0);
IndexAgentsDontLeaveIso = (CurrentSymptomStatus(:,1)==4 & CurrentSymptomStatus(:,2) > 0);

% Determine agents who restart isolation 
IndexAgentsRestart = (CurrentSymptomStatus(:,1)==5 & CurrentSymptomStatus(:,2) <= 0);
IndexAgentsDontRestart = (CurrentSymptomStatus(:,1)==5 & CurrentSymptomStatus(:,2) > 0);

% Update AC
AgentCharacteristics.Symptoms(IndexAgentsDontLeaveIso,2) = AgentCharacteristics.Symptoms(IndexAgentsDontLeaveIso,2) - parameters.dt;
AgentCharacteristics.Symptoms(IndexAgentsDontRestart,2) = AgentCharacteristics.Symptoms(IndexAgentsDontRestart,2) - parameters.dt;
AgentCharacteristics.Symptoms(IndexAgentsLeaveIso,1) = 0;
AgentCharacteristics.Symptoms(IndexAgentsLeaveIso,2) = 0;
AgentCharacteristics.Symptoms(IndexAgentsRestart,1) = 3;
AgentCharacteristics.Symptoms(IndexAgentsRestart,2) = parameters.IsolationClearanceTestDay;

         
