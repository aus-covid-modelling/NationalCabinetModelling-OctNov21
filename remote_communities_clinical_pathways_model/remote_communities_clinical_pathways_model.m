function [outcomes_by_age_and_vac,total_outcomes]=remote_communities_clinical_pathways_model(file_name,age_shift)
% Last modified: 5/11/2021 
%
% Written in: MATLAB 2019b
%
% Inputs:
% file_name - the file name of the .csv output from the agent-based
% epidemic model
% age_shift - takes values 0, 10 or 20 (depending on level of severity for
% simulation)
%
%Outputs:
% outcomes_by_age - rows correspond to number symptomatic, admitted to hospital and
% admitted to ICU. Columns refer to the cohort as <15 vac, <15 unvac, 15-39 vac, 15-39
% unvac, 40-59 vac, 40-59 unvac, 60+ vac, 60+unvac

tic;

sim_table = readtable(file_name,'ReadVariableNames',true);

assymp_symp = unique(sim_table.Symptomatic);

[~,n_cols]=size(sim_table);
sim_table=sim_table(strcmp(sim_table.Symptomatic,assymp_symp(2)),[1:3,5:n_cols]);


non_case_columns = 4;

days_full=sim_table.DateSymptoms;
day_vec=unique(days_full);
num_days = length(day_vec);

age_full = sim_table.AgeBracket;
age_strata = unique(age_full);
age_strata = {age_strata{1},age_strata{10},age_strata{[2:9,11:(end)]}}; %CHECK THIS ORDERING IS CORRECT
vacc_full = strcat(sim_table.VaccineAtInfection,num2str(sim_table.DosesAtInfection));
vacc_strata = unique(vacc_full);
num_vacc = length(vacc_strata);
num_ages = length(age_strata);
vacc_strata = vacc_strata([5,1,2,6,7,3,4]) %AlWAYS CHECK THIS ORDERING WORKS

num_strata = length(vacc_strata)*length(age_strata);

[~,x]=size(sim_table);
num_sims = x-non_case_columns;


%% Recording symptomatic cases and potential delays to presentation
cases_matrix=zeros(num_sims,num_strata,num_days);
for ii=1:num_days
    for vacc = 1:length(vacc_strata)
        for age = 1:length(age_strata)
            strata_day_index = find(strcmp(vacc_full,vacc_strata{vacc}) & strcmp(age_full,age_strata{age}) & days_full==day_vec(ii));
            cases_matrix(:,age + (vacc-1)*length(age_strata),ii) = table2array(sim_table(strata_day_index,(non_case_columns+1):end));
        end
    end
    
end
cases_matrix = permute(cases_matrix, [3,1,2]);



min_samples = 200;
increase_sample=ceil(min_samples/num_sims);
if increase_sample>1
    cases_matrix = repmat(cases_matrix,[1,increase_sample,1]);
    num_sims = increase_sample*num_sims;
end
display(['the simulations were augmented to ', num2str(num_sims)])


%% Computing outcome probability vectors 

%%Symptom Probabilities 
rel_symp_vacc = 1-[0, 0.4, 0.71, 0.58, 0.84, 0.58, 0.84];

symps_baseline = [0.28,0.28,0.2,0.2,0.26,0.26,0.33,0.33,0.4,0.4,0.49,0.49,0.63,0.63,0.69,0.69,0.69];

%computing probability of symptoms given infection, age and vaccine status
VE_temp = rel_symp_vacc;
temp = symps_baseline;
temp_p_symp = zeros(1,num_strata);
for vvv=1:num_vacc
    temp_p_symp((1:num_ages)+num_ages*(vvv-1))=min((VE_temp(vvv).*temp./(1-temp))./(1+VE_temp(vvv).*(temp./(1-temp))),1);
end
p_symp = temp_p_symp;

%%Hospitalisation Probabilities

%relative risk of VOC compared to Wuhan strain
OR = 2.08;
RR = 1;

rel_hosp_vacc = 1-[0,0.81,0.77,0.92,0.93,0.93,0.97];


if age_shift==0
    hosp_baseline = 0.75 * [0.039 0.001  0.006  0.009 0.026 0.040 0.042 0.045 0.050 0.074 0.138 0.198 0.247 0.414 0.638 1.000 0.873];
elseif age_shift==10
    hosp_baseline = 0.75 * [0.039 0.001  0.006  0.009 0.042 0.045 0.050 0.074 0.138 0.198 0.247 0.414 0.638 1.000 1.000 1.000 1.000];
elseif age_shift==20
    hosp_baseline = 0.75 * [0.039 0.001  0.006  0.009 0.050 0.074 0.138 0.198 0.247 0.414 0.638 1.000 1.000 1.000 1.000 1.000 1.000];
end

% computing probability of hospitalisation given symptoms
VE_temp = rel_hosp_vacc;
temp = hosp_baseline.*symps_baseline;
temp_p_hosp = zeros(1,num_strata);
for vvv=1:num_vacc
    temp_p_hosp((1:num_ages)+num_ages*(vvv-1)) = min((VE_temp(vvv).*temp./(1-temp))./(1+VE_temp(vvv).*(temp./(1-temp))),1);
end
p_hosp_temp = min((OR.*temp_p_hosp./(1-temp_p_hosp))./(1+OR.*(temp_p_hosp./(1-temp_p_hosp))),1);
p_hosp = min(RR * p_hosp_temp./p_symp,1);


%%ICU probabilities
voc_rel_ICU = 1; 
OR_ICU_delta=3.35;
rel_vacc_ICU = 1-[0,0.81,0.77,0.92,0.93,0.93,0.97];

% 0.24 (0.14, 0.36)
if age_shift==0
    ICU_baseline =0.24 * [0.243 0.289 0.338 0.389 0.443 0.503 0.570 0.653 0.756 0.866 0.954 1.000 0.972 0.854 0.645 0.402 0.107];
elseif age_shift==10
    ICU_baseline =0.24 * [0.243 0.289 0.338 0.389  0.570 0.653 0.756 0.866 0.954 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];
elseif age_shift==20
    ICU_baseline =0.24 * [0.243 0.289 0.338 0.389 0.756 0.866 0.954 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];
end


VE_temp = rel_vacc_ICU;
temp_p_preICU=zeros(1,num_strata);
temp=voc_rel_ICU * ICU_baseline .*hosp_baseline.*symps_baseline;
for vvv=1:num_vacc
    temp_p_preICU((1:num_ages)+num_ages*(vvv-1))=min((VE_temp(vvv).*temp./(1-temp))./(1+VE_temp(vvv).*(temp./(1-temp))),1);
end
temp2 = min((OR_ICU_delta.*temp_p_preICU./(1-temp_p_preICU))./(1+OR_ICU_delta.*(temp_p_preICU./(1-temp_p_preICU))),1);
p_preICU = min(temp2./(p_hosp.*p_symp),1);


%% Sampling Outcomes
fes_cases = squeeze(sum(cases_matrix,1));
fes_hosp=binornd(fes_cases,repmat(p_hosp,[num_sims,1]));
fes_icu=binornd(fes_hosp,repmat(p_preICU,[num_sims,1]));

%% Formatting Outcomes
fes_unvacc_cases = [sum(fes_cases(:,1:3),2),sum(fes_cases(:,4:8),2),sum(fes_cases(:,9:12),2),sum(fes_cases(:,13:num_ages),2)];
fes_unvacc_hosp = [sum(fes_hosp(:,1:3),2),sum(fes_hosp(:,4:8),2),sum(fes_hosp(:,9:12),2),sum(fes_hosp(:,13:num_ages),2)];
fes_unvacc_icu = [sum(fes_icu(:,1:3),2),sum(fes_icu(:,4:8),2),sum(fes_icu(:,9:12),2),sum(fes_icu(:,13:num_ages),2)];

fes_vacc_cases = zeros(size(fes_unvacc_cases));
fes_vacc_hosp = zeros(size(fes_unvacc_hosp));
fes_vacc_icu = zeros(size(fes_unvacc_icu));
for ii=2:num_vacc
    shift = num_ages*(ii-1);
    fes_vacc_cases = fes_vacc_cases + [sum(fes_cases(:,(1:3)+shift),2),sum(fes_cases(:,(4:8)+shift),2),sum(fes_cases(:,(9:12)+shift),2),sum(fes_cases(:,(13:num_ages)+shift),2)];
    fes_vacc_hosp = fes_vacc_hosp + [sum(fes_hosp(:,(1:3)+shift),2),sum(fes_hosp(:,(4:8)+shift),2),sum(fes_hosp(:,(9:12)+shift),2),sum(fes_hosp(:,(13:num_ages)+shift),2)];
    fes_vacc_icu = fes_vacc_icu + [sum(fes_icu(:,(1:3)+shift),2),sum(fes_icu(:,(4:8)+shift),2),sum(fes_icu(:,(9:12)+shift),2),sum(fes_icu(:,(13:num_ages)+shift),2)];
end

outcomes_by_age_and_vac = [mean(fes_vacc_cases),mean(fes_unvacc_cases);mean(fes_vacc_hosp),mean(fes_unvacc_hosp);mean(fes_vacc_icu),mean(fes_unvacc_icu)];
outcomes_by_age_and_vac = outcomes_by_age_and_vac(:,[1,5,2,6,3,7,4,8]);
total_outcomes = [sum(mean(fes_unvacc_cases + fes_vacc_cases));sum(mean(fes_unvacc_hosp + fes_vacc_hosp));sum(mean(fes_unvacc_icu + fes_vacc_icu))];

