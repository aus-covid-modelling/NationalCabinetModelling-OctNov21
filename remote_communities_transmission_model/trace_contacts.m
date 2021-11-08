function [AgentCharacteristics, parameters] = ...
    trace_contacts(parameters,AgentCharacteristics,CurrentContacts,...
    CurrentSymptomStatus, CurrentInfection, AgentIndex, KeepIndexAllAgents, CurrentIDs, timestep)

% This script back traces contacts of confirmed cases

ContactsAll = cell2mat(CurrentContacts(AgentIndex,1));
NumContactsEachAgent = cellfun('size',CurrentContacts(AgentIndex,1),1);

if sum(NumContactsEachAgent)~=0
    TimesAll = ContactsAll(:,3); % units timesteps


    % Minimum time of tracing depends on symptom status
    % Asympt: time of test - x days, Sympt: time of symptom onset - x days
    InfectionStatusIsolated = CurrentInfection(AgentIndex,1);
    TSSO = zeros(size(InfectionStatusIsolated)) + parameters.DelayOnsetSymptTestResult;  % units days
    TSSO(InfectionStatusIsolated==3) = CurrentSymptomStatus(InfectionStatusIsolated==3,3); % units days
    MinimumTimeEachAgent = timestep - (TSSO(:) + parameters.LengthofTracing) / parameters.dt; % units timesteps

    % Need to compare times of each contact to MinimumTimeEachAgent, so form
    % array where MinimumTimeEachAgent is repeated for every contact of each
    % agent

    MinimumTime = repelem(MinimumTimeEachAgent, NumContactsEachAgent);
 
    TraceableContacts = (TimesAll(:) > MinimumTime(:));
    TraceableContacts = ContactsAll(TraceableContacts,1);
    if ~isempty(TraceableContacts)
        %TraceableContacts = unique(TraceableContacts); %SLOW
        [temp, ~] = sort(TraceableContacts); %FAST
        TraceableContacts = temp([true;diff(temp(:))>0]);
    else
        TraceableContacts=[];
    end

    RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+length(TraceableContacts)-1);
    TraceableContacts = TraceableContacts(RandomNumbers<parameters.EffectivenessContactTracing);

    % Update counter for pregenerated random numbers:
    parameters.countUR = parameters.countUR + length(TraceableContacts);
    if parameters.countUR > 0.9*10^5
        parameters.UniformRand = rand(1e5,1);
        parameters.countUR = 1;
    end

    % Remove agents who are in isolation (for those in self isolation, keep in quarantine pathway)
    Isolated = KeepIndexAllAgents(CurrentSymptomStatus(:,1)>1);
    IsolatedIDs = CurrentIDs(Isolated);
    TraceableContacts(ismember(TraceableContacts,IsolatedIDs))=[];

    % Remove agents who are not currently in population
    TraceableContacts(~ismember(TraceableContacts,CurrentIDs))=[];

    % Find row indices of remaining contacts 
    [~,RowIndexesTraceableContacts]=ismember(TraceableContacts,CurrentIDs);

    % Start contacts in quarantine pathway
    AgentCharacteristics.Quarantine(RowIndexesTraceableContacts,1) = 1;
    AgentCharacteristics.Quarantine(RowIndexesTraceableContacts,2) = parameters.DelayCaseIdentifiedQuarantine;
end
