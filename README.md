# gazeHMM Validation
This is the online repository for the paper "Classifying Eye Movement Events With an Unsupervised Generative Hidden Markov Model". A preprint version of the paper can be found here: https://doi.org/10.31234/osf.io/wvp2f. This repository contains the online supplementary material for the paper as well as the code to reproduce the results and the manuscript.

## Structure

- *algorithm*: R functions for gazeHMM algorithm  
- *manuscript*: files for compiling the preprint manuscript, supplementary material
- *simulation*: R scripts for running the simulation study  
  + *preregistration*: files for compiling the preregistration of the simulation study  
- *validation*: R scripts for running the validation of gazeHMM  

## Reproduction
The preprint manuscript can be reproduced by running *preprint_Luken_Kucharsky_Visser_Classifying_Eye_Movement_Events.Rmd*. Several files are required for the reproduction:

- Simulation results - can be obtained by running *parameter_recovery_simulation.R* and *parameter_recovery_simulation_exploration.R*; the simulation takes a lot of time to run and thus, the results are included in the repository, i.e., *part_X.Rdata* and *part_3_expl.Rdata*; an image of R after the simulation was run is contained in *results_image.Rdata*
- Raw data and fitted algorithm data for the Andersson et al. (2017) data set: Those can be obtained by placing the data of the original article in *validation/data* and running *validation_Andersson2017.R*
- Fitted algorithm data for the Ehinger et al. (2019) data set, which can be obtained by placing the .EDF (EyeLink) files of the original article in *validation/data* and running *validation_Ehinger2019.R*

## References
The references for the two validation data sets are:

Andersson, R., Larsson, L., Holmqvist, K., Stridh, M., & Nyström, M. (2017). One algorithm to rule them all? An evaluation and discussion of ten eye movement event-detection algorithms. *Behavior Research Methods, 49*, 616-637. <https://doi.org/10.3758/s13428-016-0738-9>

Ehinger, B. V., Groß, K., Ibs, I., König, P. (2019). A new comprehensive eye-tracking test battery concurrently evaluating the Pupil Labs glasses and the EyeLink 1000. *PeerJ 7*, e7086. <https://doi.org/10.7717/peerj.7086>

***
Check out the R package for gazeHMM under https://github.com/maltelueken/gazeHMM!
