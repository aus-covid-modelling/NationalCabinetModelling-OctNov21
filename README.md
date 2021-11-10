# Modelling undertaken for National Cabinet in October/November 2021

# Work Package 1
The various models for assessment of test, trace, isolation and quarantine (TTIQ) are contained in the external repository [ttiq-simulation](https://github.com/njtierney/ttiq-simulation/tree/7b9897bc37c3c577d19bdf633249210e67ad742d).

# Work Package 2

## Local Government and small area effects substream
The script to reproduce the work from LGAs is contained in [contact_matrices_australia/lgas.R](https://github.com/goldingn/contact_matrices_australia/blob/237bda4a083721801b91d63a197e86f71b98feb2/lgas.R). The work-from-home effect can be calculated using [contact_matrices_australia/wfh.R](https://github.com/goldingn/contact_matrices_australia/blob/237bda4a083721801b91d63a197e86f71b98feb2/wfh.R). Maps can be drawn using [contact_matrices_australia/mapping.R](https://github.com/goldingn/contact_matrices_australia/blob/237bda4a083721801b91d63a197e86f71b98feb2/mapping.R).

## First Nations Austalians substream
The script to produce the effects on transmission potential by area is in [contact_matrices_australia/nt_aboriginal_communities.R](https://github.com/goldingn/contact_matrices_australia/blob/237bda4a083721801b91d63a197e86f71b98feb2/nt_aboriginal_communities.R)

The transmission model and clinical pathways model are located in [remote_communities_transmission_model](remote_communities_transmission_model) and [remote_communities_clinical_pathways_model](remote_communities_clinical_pathways_model) respectively. Please see the separate READMEs for each model (located [here](remote_communities_transmission_model/README_remote_communities_abm.md) and [here](remote_communities_clinical_pathways_model_model/README_remote_communities_cp.md)) for more specific information on these models.

## Schools substream
The COVASIM model used for the Schools substream is in [schools_abm](schools_abm). Please see the readme in that repository for details on exactly how to run that model.

# Work Package 3
The quarantine model is contained in [quarantine_abm](quarantine_abm). This agent based model is what is used to generate the breach events from the quarantine system to be seeded in the community.

The population level model is contained in [population_abm_package3](population_abm_package3). This model takes outputs from the quarantine ABM, transforms them into an appropriate format (using [create_arrivals_linelist.R](population_abm_package3/create_arrivals_linelist.R)) and seeds infections in the community.

# Parameter updates
The script to update estimates of infectiousness and susceptibility by age is in [contact_matrices_australia/calibrate_settings.R](https://github.com/goldingn/contact_matrices_australia/blob/237bda4a083721801b91d63a197e86f71b98feb2/calibrate_settings.R)

The software to extrapolate age-structured contact matrices based on demographic information is available ain the repo [njtierney/conmat](https://github.com/njtierney/conmat). [This commit](https://github.com/njtierney/conmat/tree/08bfcd22b39f449f330c27e7a59f7009854250de) gives the version used in the final reports.
