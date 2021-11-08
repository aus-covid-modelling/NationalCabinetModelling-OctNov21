function [AgentCharacteristics, parameters] = initialise_demographics(parameters)

% This script initialises age, dwelling set, current dwelling, unique agent ID for all agents 

% Specify core houses of each agent
RelevantHHIDs = parameters.HHIDs;
CoreHH = datasample(1:parameters.NumberHouses,parameters.PopSize);
CoreHH = RelevantHHIDs(CoreHH);
hc = histogram(CoreHH,parameters.NumberHouses);
vc = hc.Values; % this is array of core hh sizes
close all

% Make sure no houses are without core residents.  If there are,
% add additional members to the population to populate
if isempty(vc(vc==0))==0
    NumberNewPeople = length(vc(vc==0));
    x = find(vc==0);
    CoreHH = [CoreHH; x(:)];
    parameters.PopSize = parameters.PopSize + NumberNewPeople;
    parameters.PopSize = sum(parameters.PopSize);
end

% Specify regular and OnOff houses of each agent
RegularHH = datasample(1:parameters.NumberHouses,parameters.PopSize);
RegularHH = RelevantHHIDs(RegularHH);
OnOffHH = datasample(1:parameters.NumberHouses,parameters.PopSize);
OnOffHH = RelevantHHIDs(OnOffHH);
HouseholdList = [CoreHH(:) RegularHH(:) OnOffHH(:)];

% Specify age for each agent
load('data/age_distribution.mat','B')
Age = datasample(0:1:89,parameters.PopSize,'Weights',B);
Age(Age>=parameters.AgeDeath)=0;
for i = 0:parameters.AgeDeath-1
    ind = (Age==i);
    newages = rand(sum(ind),1)+i;
    Age(ind)=newages;
end

if parameters.exemplar == 1 
    age_dist = [45 17 89 49 20 0]; 
    age_groups = [0 12 16 40 60 80];
    widths = diff([age_groups parameters.AgeDeath]);
    Age = datasample(age_groups,parameters.PopSize,'Weights',age_dist);
    for i = 1:length(age_dist)
        lowage = age_groups(i);
        ind = (Age==lowage);
        newages = rand(sum(ind),1)*widths(i)+lowage;
        Age(ind)=newages;
    end
    Age(Age>=parameters.AgeDeath)=0;
elseif parameters.exemplar == 2 
    age_dist = [155 40 234 109 39 3]; 
    age_groups = [0 12 16 40 60 80];
    widths = diff([age_groups parameters.AgeDeath]);
    Age = datasample(age_groups,parameters.PopSize,'Weights',age_dist);
    for i = 1:length(age_dist)
        lowage = age_groups(i);
        ind = (Age==lowage);
        newages = rand(sum(ind),1)*widths(i)+lowage;
        Age(ind)=newages;
    end
    Age(Age>=parameters.AgeDeath)=0;
elseif parameters.exemplar == 3 
    age_dist = [93 161 75 83 75 118 59 65 46 66 50 42 33 20 21 9 2]; 
    age_groups = [0 5 12 16 20:5:80];
    widths = diff([age_groups parameters.AgeDeath]);
    Age = datasample(age_groups,parameters.PopSize,'Weights',age_dist);
    for i = 1:length(age_dist)
        lowage = age_groups(i);
        ind = (Age==lowage);
        newages = rand(sum(ind),1)*widths(i)+lowage;
        Age(ind)=newages;
    end
    Age(Age>=parameters.AgeDeath)=0;
end

% Number of agents in each household
CurrentHousehold = HouseholdList(:,1); % start in core hh
parameters.nahh = histcounts(CurrentHousehold,1:parameters.NumberHouses+1);

ID=3+(1:1:parameters.PopSize)';
parameters.countID = max(ID);

% AgentCharacteristics.Age(i): is the age of agent i  
AgentCharacteristics.Age = Age';

% AgentCharacteristics.HouseholdList(i,:) is the [core regular on/off] hh 
% IDs of agent i
AgentCharacteristics.HouseholdList = HouseholdList;

% AgentCharacteristics.CurrentHousehold(i) is the hh ID of the house the
% agent is currently residing in
AgentCharacteristics.CurrentHousehold = CurrentHousehold;

% AgentCharacteristics.ID(i) is the unique ID of agent i
AgentCharacteristics.ID = ID;

end
