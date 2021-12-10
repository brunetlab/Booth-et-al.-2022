Here I am analyzing the lifespan results using Cox Proportional Hazard Regression. 

The scripts for this set of data analysis are:

1. Cox-Proportional-Hazard-Regression.R
This example script calculates the hazard ratios, confidence intervals, and p-values for an experiment in the manuscript. It determines the hazard ratios for the interaction between male-induced demise and a given mutation/RNAi and the hazard ratios for a given mutation/RNAi in both the her only and her (+males) conditions. The same script was used for each lifespan with slight modification for the different conditions.

2. balloon_plot_cox_screen_longevity.R
To create Fig. 2a, this script generates a balloon plot summarizing the results of the RNAi screen as well as a representative lifespan experiment testing the susceptibility of classical longevity genes to male-induced demise. The output was modified in illustrator to change the order and to omit the interaction results for simplicity.

3. box-plots.R
This script uses ggplot2 to create summary box plots (shown in Fig. 2m and Ext. Data Fig. 2a).



