function parameters = parameters_infection(parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all pathogen parameters for baseline simluation here %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disease states:
% S, E -> PI, I -> R
%      -> PA, A -> R

% Vaccine states:
% VP1, VP2, VAZ1, VAZ2

% Initialisation
% Initial number in infection states
parameters.NE0 = 1; 
parameters.NPI0 = 0;  
parameters.NI0 = 0; 
parameters.NPA0 = 0; 
parameters.NA0 = 0; 
parameters.NR0 = 0; 

% Initial proportion vaccinated (for each age groups)
% Age groups: <12,	12-<15,	15-<40,	40-<60,	60+
parameters.VP10 = [0 0 0 0 0]; % 1 dose Pfizer
parameters.VP20 = [0 0 0 0 0]; % 2 dose Pfizer
parameters.VAZ10 = [0 0 0 0 0]; % 1 dose AZ
parameters.VAZ20 = [0 0 0 0 0]; % 2 dose AZ

% Within-host parameters
parameters.Ptransmission = 0.543; % Calibrate to achieve desired R0 

% Pr(symptoms|infection) and relative susceptibility
% (Mean) Davies et al. Nature 2021
parameters.SymptomsACD = [0:10:70 parameters.AgeDeath];
parameters.ProbabilitySymptoms = [0.29
0.21
0.27
0.33
0.4
0.49
0.63
0.69];

% Relative transmission probabilities
parameters.TPACD = [0:5:80 parameters.AgeDeath];
filename = 'data/setting_transmission_probabilities_for_bec.xlsx';
sheet = 1;
xlRange = 'D2:D579';
tp = xlsread(filename,sheet,xlRange);

for i = 1:2
    tpm=zeros(17,17);
    for j = 1:17
        rows = (i-1)*17*17 + (j-1)*17+(1:17);
        tpm(j,:) = tp(rows)';
    end
    if i == 1
        parameters.RelativeTransmissibilityHH = tpm;
    else
        parameters.RelativeTransmissibilityCom = tpm;
    end
end

parameters.RelativeTransmissibilityHH = parameters.RelativeTransmissibilityHH';
parameters.RelativeTransmissibilityCom = parameters.RelativeTransmissibilityCom';

% Relative reduction in onward transmission from asymptomatic
parameters.ReductionOnwardTAsymptomatic = 0.5;

% Durations
% Latency
parameters.MeanDurationLatency = 1.3; %lognormal mu
parameters.SDLatency = 0.4; %lognormal sigma
% Incubation
parameters.MeanDurationIncubation = 1.51;  %lognormal mu
parameters.SDIncubation = 0.46; %lognormal sigma
% Symptoms (from calibration)
parameters.MeanDuratioSymptoms = 1.4010; %lognormal mu
parameters.SDSymptoms = 0.0987; %lognormal sigma

% Vaccination parameters
parameters.Max_VE_Reduction_Infection_Pfizer_Dose1 = 0.57;
parameters.Max_VE_Reduction_Infection_Pfizer_Dose2 = 0.80;
parameters.Max_VE_Reduction_SymptomaticInfection_Pfizer_Dose1 = 0.58;
parameters.Max_VE_Reduction_SymptomaticInfection_Pfizer_Dose2 = 0.84;
parameters.Max_VE_Reduction_OnwardsT_Pfizer_Dose1 = 0.13;
parameters.Max_VE_Reduction_OnwardsT_Pfizer_Dose2 = 0.65;
parameters.Max_VE_Reduction_Infection_AZ_Dose1 = 0.47;
parameters.Max_VE_Reduction_Infection_AZ_Dose2 = 0.67;
parameters.Max_VE_Reduction_SymptomaticInfection_AZ_Dose1 = 0.40;
parameters.Max_VE_Reduction_SymptomaticInfection_AZ_Dose2 = 0.71;
parameters.Max_VE_Reduction_OnwardsT_AZ_Dose1 = 0.02;
parameters.Max_VE_Reduction_OnwardsT_AZ_Dose2 = 0.36;

% Onward transmission correction factors
filename = 'data/onwards_correction_factor.xlsx';
sheet = 1;
xlRange = 'J2:J69';
tpc = xlsread(filename,sheet,xlRange);

tpcf=zeros(17,4);
for j = 1:17
    rows = (j-1)*4 + (1:4);
    tpcf(j,:) = tpc(rows)';
end
parameters.OnwardsTCorrection = tpcf;

% Dynamic vaccine efficacy for recently vaccinated
% The effects below are between [0,1], where 1 represents Max_VE_Reduction
dose1_first_wk = (0:1:(7 / parameters.dt))/(7 / parameters.dt)*0.01;
dose1_next_2wks = 0.01 + 0.99 * (1:1:(14 / parameters.dt))/(14 / parameters.dt);
dose2_2wks = (0:1:(14 / parameters.dt))/(14 / parameters.dt);

parameters.VE_Effect_Over_Time_Dose1 = [dose1_first_wk dose1_next_2wks];
parameters.VE_Effect_Over_Time_Dose2 = dose2_2wks;

% ALL OR NOTHING (not tested for reactive scenarios, only pre-emptive)
parameters.allornothingve = 0;
parameters.Max_VE_Reduction_Infection_Pfizer_Dose1_AON = 0.46;
parameters.Max_VE_Reduction_Infection_Pfizer_Dose2_AON = 0.79;
parameters.Max_VE_Reduction_Infection_AZ_Dose1_AON = 0.63;
parameters.Max_VE_Reduction_Infection_AZ_Dose2_AON = 0.93;

% Adjust second dose AON probabilities, as they will only apply to agents
% who did not reach full protection after dose 1
parameters.Max_VE_Reduction_Infection_Pfizer_Dose2_AON = (parameters.Max_VE_Reduction_Infection_Pfizer_Dose2_AON - ...
    parameters.Max_VE_Reduction_Infection_Pfizer_Dose1_AON)/(1 - ...
    parameters.Max_VE_Reduction_Infection_Pfizer_Dose1_AON);
parameters.Max_VE_Reduction_Infection_AZ_Dose2_AON = (parameters.Max_VE_Reduction_Infection_AZ_Dose2_AON - ...
    parameters.Max_VE_Reduction_Infection_AZ_Dose1_AON)/(1 - ...
    parameters.Max_VE_Reduction_Infection_AZ_Dose1_AON);

% Calculate probability per timestep (during time dependent VE stage) of 
% obtaining full protection for recentlyn vaccinated

tsd1 = length(dose1_next_2wks);
tsd2 = length(parameters.VE_Effect_Over_Time_Dose2) - 1;

parameters.ProbabilityPerTimeStep0to1VE_Pfizer_Dose1 = ...
    1 - (1 - parameters.Max_VE_Reduction_Infection_Pfizer_Dose1_AON)^(1/tsd1);
parameters.ProbabilityPerTimeStep0to1VE_Pfizer_Dose2 = ...
    1 - (1 - parameters.Max_VE_Reduction_Infection_Pfizer_Dose2_AON)^(1/tsd2);
parameters.ProbabilityPerTimeStep0to1VE_AZ_Dose1 = ...
    1 - (1 - parameters.Max_VE_Reduction_Infection_AZ_Dose1_AON)^(1/tsd1);
parameters.ProbabilityPerTimeStep0to1VE_AZ_Dose2 = ...
    1 - (1 - parameters.Max_VE_Reduction_Infection_AZ_Dose2_AON)^(1/tsd2);

