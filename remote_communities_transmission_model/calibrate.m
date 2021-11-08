function Sumstats = calibrate(mu,sigma,beta)

% This script can be used to calibrate the model.  Free parameters are
% related to the distribution for the duration of symptoms, mu and sigma,
% and the scale parameter for the probability of transmission, beta.
% Outputs a set of summary statistics:
% [R0; secondary household attack rate (note this only defined when mobility 
% bewteen dwellings is switched off); generation interval; 
% then parameters describing the distrubution of time of onset of symptoms 
% to transmission (TOST):
% Pr(TOST < -5 days); Pr(TOST < -1 days); Pr(TOST < 0 days); 
% Pr(TOST < 1 days); Pr(TOST < 5 days)]

theta = [mu,sigma,beta];

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(4); %insert max number of cores here
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GI/TOST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endtime = 0.2; % Maximum duration of each simulation in years
n = 500;

% This is where we will store the initialised agent characteristics and 
% parameters to run each simulation:
ACmaster = cell(1,n);
pmaster = cell(1,n);

t1 = theta(1);
t2 = theta(2);
t3 = theta(3);

parfor i = 1:n
    
    parameters=[];
    
    parameters.exemplar = 1;
    
    % Load population parameters
    parameters = parameters_base(parameters);

    % Load infection-related parameters
    parameters = parameters_infection(parameters);
    
    parameters.VP10 = [0 0 0 0 0];
    parameters.VP20 = [0 0 0 0 0];
    parameters.VAZ10 = [0 0 0 0 0];
    parameters.VAZ20 = [0 0 0 0 0];
        
    % Load response-related parameters
    parameters = parameters_response(parameters);
    
    parameters.response = 0;

    % Initialise age structure and household structure
    [AgentCharacteristics, parameters] = initialise_demographics(parameters);
    
    parameters.T = endtime;
    parameters.time = 0 : parameters.dt : round(parameters.T * 365 / parameters.dt);
    parameters.Ntimesteps = length(parameters.time);
    % Uncomment below two lines to calibrate model to data relevant to 
    % settings where it can be assumed there is no moblity bewteen
    % dwellings
	%parameters.HHMobilityPD = [1 0 0 0];   
	%parameters.HHMobilityPD = cumsum(parameters.HHMobilityPD);
    
    parameters.MeanDuratioSymptoms = t1;
    parameters.SDSymptoms = t2;
    parameters.Ptransmission = t3;

    % Initialise infection and immunity status
    AgentCharacteristics = initialise_infection_status(AgentCharacteristics,parameters);
    
    ACmaster{i} = AgentCharacteristics;
    pmaster{i} = parameters;

end

GT = cell(1,n);
TOST = cell(1,n);
R0 = zeros(1,n);
SARhh = zeros(1,n);

parfor j = 1:n
            
    % Load paramters and agent characteristics
    parameters = pmaster{j};
    AgentCharacteristics = ACmaster{j};

    % Simulate outbreak
    [~, SummaryStatistics, parameters] = ...
        simulator(AgentCharacteristics, parameters);
    
    % Get R0 data
    if ~isempty(SummaryStatistics.InfectorInfecteePairs)
        InitialInfecteeID = SummaryStatistics.InfectorInfecteePairs(1,1);
        NumberSusceptibleHH = SummaryStatistics.InfectorInfecteePairs(1,7);
        if ~isempty(InitialInfecteeID)
            ind = find(SummaryStatistics.InfectorInfecteePairs(:,1)==InitialInfecteeID);
            indh = find(SummaryStatistics.InfectorInfecteePairs(:,1)==InitialInfecteeID &...
                SummaryStatistics.InfectorInfecteePairs(:,4)==1);
            R00 = length(ind);
            R00h = length(indh);
            if NumberSusceptibleHH>0
                SARhh0 = R00h / NumberSusceptibleHH;
            else
                SARhh0 = 0;
            end
        else
            R00 = 0;
            SARhh0 = 0;
        end

        R0(1,j) = R00;
        SARhh(1,j) = SARhh0;
    else
        R0(1,j) = 0;
        SARhh(1,j) = 0;
    end
    
    % Get generation time data 
    if ~isempty(SummaryStatistics.InfectorInfecteePairs)
        IDsinfectors = SummaryStatistics.InfectorInfecteePairs(:,1);
        IDsinfecteds = SummaryStatistics.InfectorInfecteePairs(:,2);
        TimeInfected = SummaryStatistics.InfectorInfecteePairs(:,6);
        GenTime=[];
        for i = 1:length(IDsinfecteds)
            a = find(IDsinfectors==IDsinfecteds(i));
            if ~isempty(a)
                a = a(1);
                GenTime = [GenTime; (TimeInfected(a) - TimeInfected(i))];
                break
            end
        end

        GT{1,j} = GenTime * parameters.dt;
    
        % TOST data (find first infection for each infector)
        [IDsinfectors, ia, ~] = unique(SummaryStatistics.InfectorInfecteePairs(:,1));
        TimeTheyInfectedAnother =  SummaryStatistics.InfectorInfecteePairs(ia,6);
        TSO = SummaryStatistics.TimeSymptomOnset;

        TOSTTime=[];
        for i = 1:length(IDsinfectors)
            a = find(TSO(:,1)==IDsinfectors(i));
            if ~isempty(a)
                TOSTTime = [TOSTTime; (TimeTheyInfectedAnother(i) - TSO(a,2))];
                break
            end
        end

        TOST{1,j} = TOSTTime * parameters.dt;
    else
        GT{1,j} = {};
        TOST{1,j} = {};
    end
    
end

GenTime = [];
TOSTTime2 = [];
parfor j = 1:n
    if ~isempty(GT{1,j})
        GenTime = [GenTime;GT{1,j}]; 
    end
    if ~isempty(TOST{1,j})
        TOSTTime2 = [TOSTTime2; TOST{1,j}];
    end
end

if ~isempty(GenTime)
    GImean = mean(GenTime);
    TOSTn5 = length(TOSTTime2(TOSTTime2<-5))/length(TOSTTime2);
    TOSTn1 = length(TOSTTime2(TOSTTime2<-1))/length(TOSTTime2);
    TOST0 = length(TOSTTime2(TOSTTime2<0))/length(TOSTTime2);
    TOSTp1 = length(TOSTTime2(TOSTTime2<1))/length(TOSTTime2);
    TOSTp5 = length(TOSTTime2(TOSTTime2<5))/length(TOSTTime2);
else
    GImean = 0;
    TOSTn5 = 0;
    TOSTn1 = 0;
    TOST0 = 0;
    TOSTp1 = 0;
    TOSTp5 = 0;
end

R0mean = mean(R0);
SARhhmean = mean(SARhh);

Sumstats = [R0mean; SARhhmean; GImean; TOSTn5; TOSTn1; TOST0; TOSTp1; TOSTp5];

