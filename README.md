# Nowcasting German COVID-19 hospitalizations

This is the public code repository for the **SU** contribution to the <a href="https://github.com/KITmetricslab/hospitalization-nowcast-hub">German hospitalization nowcast hub</a>.

### Structure of repository

- The folder **/code/** contains the code to perform the nowcasting for a specific day in all strata. It also includes code to post-process resulting model fits into the submission format expected by the German hospitalization nowcast hub.
- The subfolder **/code/stan** contains two STAN models that implement the SU model used for Nowcasting in the Nowcast hub. The model is closely related to the Bayesian hierarchical nowcasting proposed in <a href="https://onlinelibrary.wiley.com/doi/10.1002/bimj.202000112">Günther et al. (2020)</a>. 
- The folder **/hospitalization-nowcast-hub/** is a reduced version of the <a href="https://github.com/KITmetricslab/hospitalization-nowcast-hub">Github repository</a> of the German hospitalization nowcast hub, containing data for a specific day: "2023-05-12". This folder exists to illustrate the analysis pipeline performed on a daily basis for the SU contribution.
- The folder **/results/** is used to store intermediate results of the strata-specific model fitting via HMC in STAN.



