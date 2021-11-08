% Generate data for clinical pathways model from transmission model data

reactive = 1; % for reactive scenarios, set as 1, for pre-emptive sceanrios, set as 0

if reactive == 1
    filename = 'data/reactive_scenarios.xlsx';
    sheet = 1;
    xlRange = 'A2:AD85';
elseif reactive == 0
    filename = 'data/preemptive_scenarios.xlsx';
    sheet = 1;
    xlRange = 'A2:Q25';
end

inputs = xlsread(filename,sheet,xlRange);

for i = 1:size(inputs,1)
    if reactive == 1
    	fname = sprintf( '%s%i%s', 'cpw_outputs/batch_reactive_', i,'.mat');
    elseif reactive == 0
        fname = sprintf( '%s%i%s', 'cpw_outputs/batch_preemptive_', i,'.mat');
    end
    load(fname)

    number_simulations = size(SSmaster,2);
    store = zeros(801*14*17,number_simulations);

    for j = 1:number_simulations
        cpw_data = SSmaster{1,j}.ClinicalPathwaysData;
        if ~isempty(cpw_data)
            max_day = max(cpw_data(:,1));
            Day_data = cpw_data(:,1);
            Age_group = cpw_data(:,2);
            Pfizer_data = cpw_data(:,3);
            AZ_data = cpw_data(:,4);
            if reactive == 1
                Mod1_data = cpw_data(:,5);
            	Mod2_data = cpw_data(:,6);
            end
            for k = 0:max_day
                for l = 1:17
                    number_Pfizer_Dose1 = sum(Day_data==k & Pfizer_data==1 & Age_group==l);
                    number_Pfizer_Dose2 = sum(Day_data==k & Pfizer_data==2 & Age_group==l);
                    number_AZ_Dose1 = sum(Day_data==k & AZ_data==1 & Age_group==l);
                    number_AZ_Dose2 = sum(Day_data==k & AZ_data==2 & Age_group==l);
                    if reactive == 1
                    	number_M_Dose1 = sum(Day_data==k & Mod1_data==2 & Age_group==l);
                    	number_M_Dose2 = sum(Day_data==k & Mod2_data==2 & Age_group==l);
                    else
                        number_M_Dose1 = 0;
                        number_M_Dose2 = 0;
                    end
                    number_0_Dose = sum(Day_data==k & AZ_data==0 & Pfizer_data==0 & Age_group==l);

                    data = [0; 0;
                        number_AZ_Dose1;
                        number_AZ_Dose2;
                        0; 0; 
                        number_M_Dose1; 
                        number_M_Dose2; 
                        0;
                        number_0_Dose;
                        0; 0;
                        number_Pfizer_Dose1;
                        number_Pfizer_Dose2];

                    rows = (((l-1)*14+1):(14*l)) + 17*14*k;

                    store(rows,j) = data;
                end

            end
        end
    end

    filename1 = 'data/remote_communities_template.csv';
    % Read the CSV as a table (this doesnt work in some MATLAB versions, it works in 2020a)
    t = readtable(filename1,'PreserveVariableNames',true);
    % Add a new column to the end of the table
    numOfColumn = size(t, 2);

    for  j = 1:number_simulations
        storej = array2table(store(:,j),'VariableNames',cellstr(string(j)));
        t = table(t,storej);
        t = splitvars(t);
    end

    if not(isfolder('cpw_outputs/clinical_pathways_model_inputs'))
        mkdir('cpw_outputs/clinical_pathways_model_inputs')
    end
    
    % Write to CSV file
    if reactive == 1
    	fname = sprintf('%s%i%s', 'cpw_outputs/clinical_pathways_model_inputs/remote_communities_cpw_reactive_', i,'.csv');
    else
        fname = sprintf('%s%i%s', 'cpw_outputs/clinical_pathways_model_inputs/remote_communities_cpw_preemptive_', i,'.csv');
    end
    writetable(t, fname)
end

