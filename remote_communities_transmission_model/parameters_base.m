function parameters = parameters_base(parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all population parameters for baseline simluation here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model setup parameters
parameters.dt = 0.25; % Time step duration (days)
parameters.dt_years = parameters.dt / 365; % Time step duration (years)
parameters.T = 2; % Number of years to run simulations
parameters.time = 0 : parameters.dt : round(parameters.T * 365); % vector of time steps
parameters.Ntimesteps = length(parameters.time); % Number of time steps

% Population parameters
parameters.PopSize = 1000;
parameters.MeanHHsize = 7.7;

if parameters.exemplar == 1 
    parameters.PopSize = 220;
    parameters.MeanHHsize = 6.1;  
elseif parameters.exemplar == 2 
    parameters.PopSize = 580;
    parameters.MeanHHsize = 4.8;
elseif parameters.exemplar == 3 
    parameters.PopSize = 1018;
    parameters.MeanHHsize = 3.5;
end

parameters.NumberHouses = round(parameters.PopSize/parameters.MeanHHsize);
parameters.HHIDs = (1:1:parameters.NumberHouses)';

% Probability agents stay in [Core, Regular, On/Off, Sporadic] residence
% each night. Core, Regular, On/Off are specified, sporadic is a random
% house in community
parameters.HHMobilityPD = [0.66 0.23 0.09 0.02];
parameters.HHMobilityPD = cumsum(parameters.HHMobilityPD);

% Age structure in the model
parameters.AgeDeath = 84; % Max age of agent;

parameters.AgeClassDividersContacts = [0:5:80 parameters.AgeDeath];
parameters.NumberAgeClassesContacts = length(parameters.AgeClassDividersContacts) - 1;
filename = 'data/remote_non_household_contact_matrix.xlsx';
sheet = 1;
xlRange = 'B2:R18';
parameters.Ncontacts = (xlsread(filename,sheet,xlRange)') * parameters.dt;
