function run_preemptive(theta, n, batchnumber)

%theta = [parameters.PopSize;
% parameters.ipolicy;
% parameters.qpolicy;
% parameters.EffectivenessLocC;
% parameters.VP10;
% parameters.VP20;
% parameters.VAZ10;
% parameters.VAZ20;
% R0 scenario;
% parameters.allornothingve]

ps = theta{1};
lp = theta{2};
qp = theta{3};
el = theta{4};
v1 = theta{5};
v2 = theta{6};
v3 = theta{7};
v4 = theta{8};
r0 = theta{9};
aon = theta{10};

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(16); %insert max number of cores here
end

% This is where we will store the initialised agent characteristics and 
% parameters to run each simulation:
ACmaster = cell(1,n);
pmaster = cell(1,n);

parfor i = 1:n
    
    parameters=[];
    
    parameters.exemplar = 0;

    % Load population parameters
    parameters = parameters_base(parameters);
    
    parameters.PopSize = ps;
    parameters.NumberHouses = round(ps/parameters.MeanHHsize);
    parameters.HHIDs = (1:1:parameters.NumberHouses)';

    % Load infection-related parameters
    parameters = parameters_infection(parameters);
    
    parameters.VP10 = v1;
    parameters.VP20 = v2;
    parameters.VAZ10 = v3;
    parameters.VAZ20 = v4;
    
    parameters.allornothingve = aon;
    
    if r0 == 10.7
        parameters.Ptransmission = 0.543;
    else
        parameters.Ptransmission = 0.242;
    end
        
    % Load response-related parameters
    parameters = parameters_response(parameters);
    
    parameters.EffectivenessLocC = el;
    parameters.ipolicy = lp;
    parameters.qpolicy = qp;
    if parameters.qpolicy == 1
        parameters.EffectivenessQuaH = 1; 
        parameters.EffectivenessQuaC = 1; 
    elseif parameters.qpolicy == 2
        parameters.EffectivenessQuaH = 0;
        parameters.EffectivenessQuaC = 1; 
    end
    
    parameters.reactivevacc = 0;

    % Initialise age structure and household structure
    [AgentCharacteristics, parameters] = initialise_demographics(parameters);
    
    % Initialise infection and immunity status
    AgentCharacteristics = initialise_infection_status(AgentCharacteristics,parameters);
    
    ACmaster{i} = AgentCharacteristics;
    pmaster{i} = parameters;

end

SSmaster = cell(1,n);

parfor j = 1:n
            
    % Load paramters and agent characteristics
    parameters = pmaster{j};
    AgentCharacteristics = ACmaster{j};

    % Simulate outbreak
    [~, SummaryStatistics, ~] = ...
        simulator(AgentCharacteristics, parameters);
    
    SSmaster{j} = SummaryStatistics
    
end

if not(isfolder('cpw_outputs'))
    mkdir('cpw_outputs')
end

fname = sprintf ( '%s%i', 'cpw_outputs/batch_preemptive_', batchnumber);
save(fname,'SSmaster')


