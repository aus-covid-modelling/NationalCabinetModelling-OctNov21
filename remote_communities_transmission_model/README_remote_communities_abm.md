{\rtf1\ansi\ansicpg1252\cocoartf2580
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 ArialMT;\f1\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Transmission model for remote Aboriginal and Torres Straight Islander communities\
\
# Running\
\
The code was written in MATLAB 2020a.\
\
\pard\pardeftab720\partightenfactor0
\cf2 \expnd0\expndtw0\kerning0
The m-file preemptive_scenarios_batch.m runs a batch of transmission and response scenarios with various preemptive vaccination set ups. Parameters defining the scenarios are in 'data/preemptive_scenarios.xlsx'. Each scenario is run n=100 times. Outputs of this script can be used to \
- create prevalence figures and data tables using generate_tables_figures_preemptive.m  \
- generate csv files for the clinical pathways model using create_cpw_data.m \
\
\pard\pardeftab720\partightenfactor0
\cf2 The m-file reactive_scenarios_batch.m runs a batch of transmission and response scenarios with various preemptive vaccination set ups. Parameters defining the scenarios are in 'data/reactive_scenarios.xlsx'. Each scenario is run n=100 times. Outputs of this script can be used to \
- create prevalence figures and data tables using generate_tables_figures_reactive.m  \
- generate csv files for the clinical pathways model using create_cpw_data.m \cf0 \kerning1\expnd0\expndtw0 \
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f1 \cf0 \
Further details of the model are provided in the Word document: Model_details.docx}