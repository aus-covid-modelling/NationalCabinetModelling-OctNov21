% This script runs a batch of scenarios with various reactive vaccination
% set ups. Parameters defining the scenarios are in 'data/reactive_scenarios.xlsx'
% Each scenario is run n times.

filename = 'data/reactive_scenarios.xlsx';
sheet = 1;
xlRange = 'A2:AD85';

inputs = xlsread(filename,sheet,xlRange);

n = 100; % Number of simulations per scenario

scenarios = 1:84;

for j = 1:size(scenarios,2)
    i = scenarios(j);
    theta{1,1} = inputs(i,2); % exemplar
    theta{1,2} = inputs(i,3); % ld policy
    theta{1,3} = inputs(i,4); % q policy
    theta{1,4} = inputs(i,5); % r policy
    theta{1,5} = inputs(i,6); % r delay
    theta{1,6} = inputs(i,7); % r vacc rate
    theta{1,7} = inputs(i,8); % ld effect
    if inputs(i,2) <=2
        theta{1,8} =  [inputs(i,9:12) 0]; % Pfizer dose 1 initial fractions
        theta{1,9} =  [inputs(i,14:17) 0]; % Pfizer dose 2 initial fractions
        theta{1,10} = [0 0 0 0 inputs(i,13)]; % AZ dose 1 initial fractions
        theta{1,11} = [0 0 0 0 inputs(i,18)]; % AZ dose 2 initial fractions
    else
        theta{1,8} =  [inputs(i,9:13) inputs(i,21:22) 0 0 0]; % Pfizer dose 1 initial fractions
        theta{1,9} =  [inputs(i,14:18) inputs(i,26:27) 0 0 0]; % Pfizer dose 2 initial fractions
        theta{1,10} = [0 0 0 0 0 0 0 inputs(i,23:25)]; % AZ dose 1 initial fractions
        theta{1,11} = [0 0 0 0 0 0 0 inputs(i,28:30)]; % AZ dose 2 initial fractions
    end
        
    theta{1,12} = inputs(i,19); %R0
    theta{1,13} = inputs(i,20); %All or nothing
    run_reactive(theta, n, inputs(i,1))   
end
