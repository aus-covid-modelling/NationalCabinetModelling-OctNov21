function [ts, parameters] = test_sensitivity(tsso, parameters) % tsso, time since symptom onset

tsso=tsso(:);

numberagents = length(tsso);
RandomNumbers = parameters.UniformRand(parameters.countUR:parameters.countUR+numberagents-1);

% Update counter for pregenerated random numbers:
parameters.countUR = parameters.countUR + numberagents;
if parameters.countUR > 0.9*10^5
    parameters.UniformRand = rand(1e5,1);
    parameters.countUR = 1;
end

C = RandomNumbers(:) .* min(tsso(:),3.5);
s = tsso(:) + C;

ts = zeros(numberagents,1);

ts(tsso <= C) = 1./(1+exp(-(1.5 + 2.2 * s(tsso <= C))));
ts(tsso > C) = 1./(1+exp(-(1.5 - 0.22 * s(tsso > C))));

    
