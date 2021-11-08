# Transmission model for remote Aboriginal and Torres Straight Islander communities

# Running

The code was written in MATLAB 2020a.

The m-file preemptive_scenarios_batch.m runs a batch of transmission and response scenarios with various preemptive vaccination set ups. Parameters defining the scenarios are in [`data/preemptive_scenarios.xlsx`](data/preemptive_scenarios.xlsx). Each scenario is run n=100 times. Outputs of this script can be used to create prevalence figures and data tables using `generate_tables_figures_preemptive.m`. Generate CSV files for the clinical pathways model using [`create_cpw_data.m`](create_cpw_data.m).

The m-file [reactive_scenarios_batch.m](reactive_scenarios_batch.m) runs a batch of transmission and response scenarios with various reactive vaccination set ups. Parameters defining the scenarios are in ['data/reactive_scenarios.xlsx'](data/reactive_scenarios.xlsx). Each scenario is run `n=100` times. Outputs of this script can be used to create prevalence figures and data tables using [generate_tables_figures_reactive.m](generate_tables_figures_reactive.m). Generate CSV files for the clinical pathways model using [`create_cpw_data.m`](create_cpw_data.m).

Further details of the model are provided in the Word document [Model_details.docx](Model_details.docx)
