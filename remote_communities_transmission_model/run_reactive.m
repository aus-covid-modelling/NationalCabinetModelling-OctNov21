function run_reactive(theta, n, batchnumber)

%theta = [parameters.exemplar;
% parameters.ipolicy;
% parameters.qpolicy;
% parameters.reactivevacc;
% parameters.DelayStartReactiveVaccination; 
% parameters.DailyMaxVaccination;
% parameters.EffectivenessLocC;
% parameters.VP10;
% parameters.VP20;
% parameters.VAZ10;
% parameters.VAZ20;
% R0 scenario;
% parameters.allornothingve];

ex = theta{1};
lp = theta{2};
qp = theta{3};
rp = theta{4};
de = theta{5};
su = theta{6};
el = theta{7};
v1 = theta{8};
v2 = theta{9};
v3 = theta{10};
v4 = theta{11};
r0 = theta{12};
aon = theta{13};


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
    
    parameters.exemplar = ex;

    % Load population parameters
    parameters = parameters_base(parameters);
    
    % Load infection-related parameters
    parameters = parameters_infection(parameters);
    
    parameters.VP10 = v1;
    parameters.VP20 = v2;
    parameters.VAZ10 = v3;
    parameters.VAZ20 = v4;
    
    parameters.allornothingve = aon;
    
    if r0 == 10.7
        if ex == 1
            parameters.Ptransmission = 0.673;%%
        elseif ex == 2
            parameters.Ptransmission = 0.7345;%%
        elseif ex == 3
            parameters.Ptransmission = 0.837; %%
        end
    else
        if ex == 1
            parameters.Ptransmission = 0.271;%%
        elseif ex == 2
            parameters.Ptransmission = 0.321;%%
        elseif ex == 3
            parameters.Ptransmission = 0.3635;%%
        end
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
    
    parameters.reactivevacc = rp;
    parameters.DelayStartReactiveVaccination = de; 
    parameters.DailyMaxVaccination = su;

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

fname = sprintf ( '%s%i', 'cpw_outputs/batch_reactive_', batchnumber);
save(fname,'SSmaster')


