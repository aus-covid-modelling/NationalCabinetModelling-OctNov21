close all
clear all

% Creates prevalence figures and saves cumulative infections in
% spreadsheet for reactive vaccination scenarios

set(0, 'DefaultFigureRenderer', 'painters');

set(0,'DefaultTextFontName','Arial')
set(0,'DefaultTextFontSize',32)
set(0,'DefaultAxesFontSize',32)
set(0,'DefaultAxesFontName','Arial')

fignames = 'cpw_figures/prevalence_reactive_';
fignames2 = 'cpw_figures/VNV_prevalence_reactive_';

filename = 'data/reactive_scenarios.xlsx';
sheet = 1;
xlRange = 'A2:AD43';

inputs = xlsread(filename,sheet,xlRange);
NumberScenarios = size(inputs,1);

% Scenarios ordered so I can copy paste from csv into report tables
% Scenarios comparing R0=10.7 vs R0=5, CTP1, exemplar 1
output_scen1 = [1 22 2 23 4 25 6 27];
left1 = [1 2 4 6]; %CTP1, for R0=10.7, for increasing reactive vacc rates
right1 = [22 23 25 27]; %CTP1, for R0=5, for increasing reactive vacc rates

% Scenarios comparing R0=10.7 vs R0=5, CTP1, exemplar 2
output_scen2 = [8 29 9 30 11 32 13 34];
left2 = [8 9 11 13]; %CTP1, for R0=10.7, for increasing reactive vacc rates
right2 = [29 30 32 34]; %CTP1, for R0=5, for increasing reactive vacc rates

% Scenarios comparing R0=10.7 vs R0=5, CTP1, exemplar 3
output_scen3 = [15 36 16 37 18 39 20 41];
left3 = [15 16 18 20]; %CTP1, for R0=10.7, for increasing reactive vacc rates
right3 = [36 37 39 41]; %CTP1, for R0=5, for increasing reactive vacc rates

output_row_v1 = 1:2:8;
output_row_nv1 = 2:2:8;

output_row_v2 = 9:2:16;
output_row_nv2 = 10:2:16;

output_row_v3 = 17:2:24;
output_row_nv3 = 18:2:24;

output_scen = [output_scen1' output_scen2' output_scen3'];
left = [left1' left2' left3'];
right = [right1' right2' right3'];
output_row_v = [output_row_v1' output_row_v2' output_row_v3'];
output_row_nv = [output_row_nv1' output_row_nv2' output_row_nv3'];

Output = zeros(24,30);
Ntimesteps = 730;

% For each comparison scenario
for k = 1:3
    % For each vaccination rate and transmission scenario
    for i = 1:8

        scenario = output_scen(i,k);
        % Load model output
        fname = sprintf ( '%s%i%s', 'cpw_outputs/batch_reactive_', scenario,'.mat');
        load(fname)
      
        if k == 1
            PopSize = 220;
        elseif k ==2
            PopSize = 580;
        else 
            PopSize = 1018;
        end
        
        NumberSimulations = size(SSmaster,2);
        
        % Initialise storage of prevalence in whole population, in
        % vaccinated subpopulation, and in unvaccinated subpopulation
        P = zeros(Ntimesteps,NumberSimulations); %prevalence of infection
        PVacc = P;
        PNonVacc = P;
        
        % Initialise storage of cumulative infections in vacc/non vacc 
        % subpopulations, for five age groups
        CI = zeros(NumberSimulations,1);
        CNV1 = CI;
        CNV2 = CI;
        CNV3 = CI;
        CNV4 = CI;
        CNV5 = CI;
        CV1 = CI;
        CV2 = CI;
        CV3 = CI;
        CV4 = CI;
        CV5 = CI;
        
        % For each simulation 
        for simulation = 1:NumberSimulations   

            if ~isempty(SSmaster{1,simulation})
                
                % Calculated whole population prevalence
                ptemp = sum(squeeze(sum(squeeze(sum(SSmaster{1,simulation}.NumberInfectionVaccinationAgeStatusTime(1:5,:,:,:),1)),1)),1);
                P(1:length(ptemp),simulation) = ptemp;
                P(:,simulation) = P(:,simulation) / PopSize * 100;
                
                % Count number of agents vaccinated/not vaccinated over
                % time
                NumberVaccinated = sum(squeeze(sum(squeeze(sum(SSmaster{1,simulation}.NumberInfectionVaccinationAgeStatusTime(:,1:4,:,:),1)),1)),1);
                NumberNotVaccinated= sum(squeeze(sum(squeeze(sum(SSmaster{1,simulation}.NumberInfectionVaccinationAgeStatusTime(:,5,:,:,:),1)),1)),1);
                
                % Use this to calculate prevalence in vaccinated and
                % unvaccinated subpopulations
                pvtemp = sum(squeeze(sum(squeeze(sum(SSmaster{1,simulation}.NumberInfectionVaccinationAgeStatusTime(1:5,1:4,:,:),1)),1)),1);
                PVacc(1:length(pvtemp),simulation) = pvtemp./ NumberVaccinated * 100;
                pnvtemp = sum(squeeze(sum(squeeze(sum(SSmaster{1,simulation}.NumberInfectionVaccinationAgeStatusTime(1:5,5,:,:),1)),1)),1);
                PNonVacc(1:length(pnvtemp),simulation) = pnvtemp./ NumberNotVaccinated * 100;

                % Count cumulative infections in vaccinated and unvaccinated subpopulations, in each age groups
                ages1 = 1:12;
                CNV1(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages1,5,:),1)),1)));
                CV1(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages1,1:4,:),1)),1)));

                ages2 = 13:15;
                CNV2(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages2,5,:),1)),1)));
                CV2(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages2,1:4,:),1)),1)));

                ages3 = 16:40;
                CNV3(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages3,5,:),1)),1)));
                CV3(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages3,1:4,:),1)),1)));

                ages4 = 41:60;
                CNV4(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages4,5,:),1)),1)));
                CV4(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages4,1:4,:),1)),1)));

                ages5 = 61:84;
                CNV5(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages5,5,:),1)),1)));
                CV5(simulation,1) = sum(squeeze(sum(squeeze(sum(SSmaster{1, simulation}.TotalInfectedEachAgeVaccGroup(ages5,1:4,:),1)),1)));

            end

        end
        
        % Now calculate quantiles for cumulative infections
        Output_CNV1 = quantile(CNV1,[0.50 0.25 0.75]);
        Output_CNV2 = quantile(CNV2,[0.50 0.25 0.75]);
        Output_CNV3 = quantile(CNV3,[0.50 0.25 0.75]);
        Output_CNV4 = quantile(CNV4,[0.50 0.25 0.75]);
        Output_CNV5 = quantile(CNV5,[0.50 0.25 0.75]);
        Output_CV1 = quantile(CV1,[0.50 0.25 0.75]);
        Output_CV2 = quantile(CV2,[0.50 0.25 0.75]);
        Output_CV3 = quantile(CV3,[0.50 0.25 0.75]);
        Output_CV4 = quantile(CV4,[0.50 0.25 0.75]);
        Output_CV5 = quantile(CV5,[0.50 0.25 0.75]);
        
        % Store in Output in relevant cells
        if ismember(scenario, left(:,k))
            [~, b] = ismember(scenario, left(:,k));
            Output((k-1)*8+2*b-1,1:15) = [Output_CV1 Output_CV2 Output_CV3 Output_CV4 Output_CV5];
            Output((k-1)*8+2*b,1:15) = [Output_CNV1 Output_CNV2 Output_CNV3 Output_CNV4 Output_CNV5];
        else
            [~, b] = ismember(scenario, right(:,k));
            Output((k-1)*8+2*b-1,16:30) = [Output_CV1 Output_CV2 Output_CV3 Output_CV4 Output_CV5];
            Output((k-1)*8+2*b,16:30) = [Output_CNV1 Output_CNV2 Output_CNV3 Output_CNV4 Output_CNV5];
        end
        
        % Now calculate quantiles for prevalence of infection over time 
        % (for figures)
        QP = quantile(P,[0.025 0.25 0.50 0.75 0.975],2)';
        y1 = QP(3,:); %median
        x1 = 1:1:Ntimesteps;
        x_axis = x1;
        x_plot =[x_axis fliplr(x_axis)];
        y_plot=[QP(2,:), flipud(QP(4,:)')']; % shading IQR

        % Prevalence whole population (plot median and IQR over time)
        figure
        plot(x_axis, y1, 'color',[0 0.5 0],'linestyle','-', 'linewidth', 2)
        hold on
        fill(x_plot, y_plot, 1,'facecolor', [0 0.5 0], 'edgecolor', 'none', 'facealpha', 0.2);
        xlim([0 150])
        ylim([0 50])
        xlabel('Time (days)')
        ylabel('Prevalence (%)')
        drawnow
        filetitle = sprintf ( '%s%i%s', fignames,scenario, '.fig');
        savefig(filetitle)
        filetitle = sprintf ( '%s%i%s', fignames,scenario);
        saveas(gcf,filetitle,'svg')
        saveas(gcf,filetitle,'pdf')

        % Prevalence vacc/non vacc (plot median and IQR over time)
        figure
        for j = 1:2

            if j==1
                % calculate quantiles of prevalence in vaccinated 
                % subpopulation over time
                QP = quantile(PVacc,[0.025 0.25 0.50 0.75 0.975],2)';
                y1 = QP(3,:); %median
                x1 = 1:1:Ntimesteps;
                x_axis = x1;
                x_plot =[x_axis fliplr(x_axis)];
                y_plot=[QP(2,:), flipud(QP(4,:)')']; % shading IQR
                plot(x_axis, y1, 'color', 'blue','linestyle','-','linewidth', 2)
                hold on
                fill(x_plot, y_plot, 1,'facecolor', 'blue' , 'edgecolor', 'none', 'facealpha', 0.2);
                hold on
            else   
                % calculate quantiles of prevalence in unvaccinated 
                % subpopulation over time
                QP = quantile(PNonVacc,[0.025 0.25 0.50 0.75 0.975],2)';
                y1 = QP(3,:); %median
                x1 = 1:1:Ntimesteps;
                x_axis = x1;
                x_plot =[x_axis fliplr(x_axis)];
                y_plot=[QP(2,:), flipud(QP(4,:)')']; % shading IQR
                plot(x_axis, y1, '-r', 'linewidth', 2)
                hold on
                fill(x_plot, y_plot, 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.2);
            end
        end
        legend('Vaccinated','','Unvaccinated','')
        xlim([0 150])
        ylim([0 50])
        xlabel('Time (days)')
        ylabel('Prevalence (%)')
        drawnow
        filetitle = sprintf ( '%s%i%s', fignames2,scenario, '.fig');
        savefig(filetitle)
        filetitle = sprintf ( '%s%i%s', fignames2,scenario);
        saveas(gcf,filetitle,'svg')
        saveas(gcf,filetitle,'pdf')

    end
end

Rownames = {'Option 1: <12 0.5', 'Option 1: <12 0.25', 'Option 1: <12 0.75',...
            'Option 1: 12-<15 0.5','Option 1: 12-<15 0.25','Option 1: 12-<15 0.75', ...
            'Option 1: 15-<40 0.5','Option 1: 15-<40 0.25','Option 1: 15-<40 0.75',...
            'Option 1: 40-<60 0.5','Option 1: 40-<60 0.25','Option 1: 40-<60 0.75',...
            'Option 1: 60+ 0.5','Option 1: 60+ 0.25','Option 1: 60+ 0.75',...
            'Option 2: <12 0.5','Option 2: <12 0.25','Option 2: <12 0.75',...
            'Option 2: 12-<15 0.5','Option 2: 12-<15 0.25','Option 2: 12-<15 0.75',...
            'Option 2: 15-<40 0.5','Option 2: 15-<40 0.25','Option 2: 15-<40 0.75',...
            'Option 2: 40-<60 0.5','Option 2: 40-<60 0.25','Option 2: 40-<60 0.75',...
            'Option 2: 60+ 0.5','Option 2: 60+ 0.25','Option 2: 60+ 0.75'};

 A = array2table(Output,'VariableNames',Rownames);
 filename = 'cpw_outputs/cumulative_infections_reactive.xlsx';
 writetable(A,filename)




