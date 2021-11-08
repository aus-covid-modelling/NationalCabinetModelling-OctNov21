% This script runs a batch of scenarios with various pre-emptive vaccination
% set ups. Parameters defining the scenarios are in 'data/preemptive_scenarios.xlsx'
% Each scenario is run n times.

filename = 'data/preemptive_scenarios.xlsx';
sheet = 1;
xlRange = 'A2:Q25';

inputs = xlsread(filename,sheet,xlRange);

n = 100; % Number of simulations per scenario

scenarios = 1:24;

for j = 1:size(scenarios,2)
    i = scenarios(j);
    theta{1,1} = inputs(i,2); % population size
    theta{1,2} = inputs(i,3); % ld policy
    theta{1,3} = inputs(i,4); % q policy
    theta{1,4} = inputs(i,5); % eff ld
    theta{1,5} =  [inputs(i,6:9) 0]; % Pfizer dose 1 initial fractions
    theta{1,6} =  [inputs(i,11:14) 0]; % Pfizer dose 2 initial fractions
    theta{1,7} = [0 0 0 0 inputs(i,10)]; % AZ dose 1 initial fractions
    theta{1,8} = [0 0 0 0 inputs(i,15)]; % AZ dose 2 initial fractions
    theta{1,9} = inputs(i,16); %R0 scenario
    theta{1,10} = inputs(i,17); %All or nothing
    run_preemptive(theta, n, inputs(i,1))   
end
