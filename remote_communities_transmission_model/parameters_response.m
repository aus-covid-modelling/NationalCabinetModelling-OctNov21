function parameters = parameters_response(parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all response parameters here                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters.response = 1; % switch response on/off

% Quarantine, isolation, lockdown, testing parameters (note: 1/2 do not 
% correspond to CTP1/2)
parameters.qpolicy = 1;
% 1: quarantined and isolated removed from community (cannot make contacts, CTP2)
% 2: isolated removed from community (quarantined can contact core hh members, CTP1)

parameters.ipolicy = 1;
% 1: lockdown has specified duration

% Durations (days)
parameters.MinDurationQua = 14; 
parameters.MinDurationIso = 10;
parameters.MinDurationLoc = 14;

% Delays (days)
parameters.DelayOnsetSymptTestResult = 1;
parameters.DelayTestResultIsolation = 1;
parameters.DelayCaseIdentifiedQuarantine = 1;
parameters.DelayCaseIdentifiedLockdown = 1;
parameters.DelayLockdownTestResultsCommunity = 2;

% Release (day)
parameters.IsolationClearanceTestDay = parameters.MinDurationIso - 2;
parameters.QuarantineClearanceTestDay = parameters.MinDurationQua - 2;
parameters.LockdownClearanceTestDay = parameters.MinDurationLoc - 2;

% Contact tracing
parameters.EffectivenessContactTracing = 1; % proportion of contacts found by tracers
parameters.LengthofTracing = 3; % Length of time can backtrace contacts from specified time

% Effectiveness Q/I/L (proportional reduction in contact rates, in hh, comm)
parameters.EffectivenessSelfIsoH = 1; % agents self isolating
parameters.EffectivenessSelfIsoC = 1; % agents self isolating
parameters.EffectivenessIsoH = 1;  % agents isolating
parameters.EffectivenessIsoC = 1;  % agents isolating

% Agents quarantining
if parameters.qpolicy == 1 % CTP2
    parameters.EffectivenessQuaH = 1; 
    parameters.EffectivenessQuaC = 1; 
elseif parameters.qpolicy == 2 % CTP1
    parameters.EffectivenessQuaH = 0;
    parameters.EffectivenessQuaC = 1; 
end

% Agents locked down
parameters.EffectivenessLocH = 0; 
parameters.EffectivenessLocC = 0.9; 

%Test sensitivity
parameters.TestSensitivity = [0 1 1 1 1 1 0];

%Lockdown switched on
parameters.lockdown = 1;

% Reactive vaccination switched on (remember to do this for "no program" 
% scenarios, and set max rate to 0)
parameters.reactivevacc = 1;
% Assume surge rate and supply are independent of schedule supply, and only 
% applies to first dose
parameters.DelayStartReactiveVaccination = 2; 
parameters.DailyMaxVaccination = 30;
parameters.Hesitancy = 0.0687;
parameters.OlderFirst = 1;

parameters.DoseIntervalP = 21;
parameters.DoseIntervalAZ = 28;

parameters.MinAgeVaccination = 12;
TimeToSeverityEffectDose1 = 14;
TimeToSeverityEffectDose2 = 5;

% Translate these times to proportional VE effect:
TimeToSD1 = TimeToSeverityEffectDose1/parameters.dt+1;
TimeToSD2 = TimeToSeverityEffectDose2/parameters.dt+1;
parameters.MinVEForSeverityDose1 = parameters.VE_Effect_Over_Time_Dose1(TimeToSD1);
parameters.MinVEForSeverityDose2 = parameters.VE_Effect_Over_Time_Dose2(TimeToSD2);