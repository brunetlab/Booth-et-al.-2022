#https://www.r-bloggers.com/examples-using-r-randomized-block-design/

rm(list=ls())
library(ggplot2)
library(survival)

##library(rms)
#library(RColorBrewer)
#https://cran.r-project.org/web/packages/qwraps2/vignettes/summary-statistics.html
#library(qwraps2)
options(stringsAsFactors = FALSE)

setwd("~/Lifespan-Analysis/")

#################################
# Organizing and Cleaning Data ##
#################################

# Loading data
df2 = read.delim("example_lifespan_data.txt", header = TRUE)

# Making sure these are factors
df2$genotype_RNAi = as.factor(df2$genotype_RNAi)
df2$males = as.factor(df2$males)


##################################
# Testing for interaction terms ##
##################################
# relevel such that baseline is control genotype_RNAi, no males
df2$genotype_RNAi = relevel(df2$genotype_RNAi, "control")
df2$males = relevel(df2$males, "no")

# Testing for genotype_RNAi*males interaction term, first with cox proportional hazards model
coxfit_df2 = coxph(Surv(Age, Status) ~ genotype_RNAi + males + genotype_RNAi:males, data = df2)
summary(coxfit_df2)


# Save results
capture.output(summary(coxfit_df2),
               file = "coxfit_results.txt")


## Create a table with the results I will report (hazard ratio, 95% confidence interval and p-value):
# Prepare the columns
HR <- exp(coef(coxfit_df2))
CI <- exp(confint(coxfit_df2))
P <- coef(summary(coxfit_df2))[,5]
# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")
# Bind columns together as dataset
results.table <- as.data.frame(cbind(HR, CI, P))
write.table(results.table,
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "coxfit_results_table.txt")



###########################
## Only +males condition ## 
###########################

MID <- df2[df2$males == "yes", ]
MID$males = as.factor(MID$males)

# relevel such that baseline is control RNAi
MID$genotype_RNAi = relevel(MID$genotype_RNAi, "control")

# Testing for effect of RNAi using cox proportional hazards model
coxfit_MID = coxph(Surv(Age, Status) ~ genotype_RNAi, data = MID)
summary(coxfit_MID)

# Save results
capture.output(summary(coxfit_MID),
               file = "coxfit_results_males-only.txt")

## Create a table with the results I will report (hazard ratio, 95% confidence interval and p-value):
# Prepare the columns
HR <- exp(coef(coxfit_MID))
CI <- exp(confint(coxfit_MID))
P <- coef(summary(coxfit_MID))[,5]
# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")
# Bind columns together as dataset
results.table <- as.data.frame(cbind(HR, CI, P))
write.table(results.table,
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "coxfit_results_table_males-only.txt")

########################
## Her only condition ##
########################

her.only <- df2[df2$males == "no", ]
her.only$males = as.factor(her.only$males)

# relevel such that baseline is control RNAi
her.only$genotype_RNAi = relevel(her.only$genotype_RNAi, "control")

# Testing for effect of RNAi using cox proportional hazards model
coxfit_her.only = coxph(Surv(Age, Status) ~ genotype_RNAi, data = her.only)
summary(coxfit_her.only)

# Save results
capture.output(summary(coxfit_her.only),
               file = "coxfit_results_her-only.txt")

## Create a table with the results I will report (hazard ratio, 95% confidence interval and p-value):
# Prepare the columns
HR <- exp(coef(coxfit_her.only))
CI <- exp(confint(coxfit_her.only))
P <- coef(summary(coxfit_her.only))[,5]
# Names the columns of CI
colnames(CI) <- c("Lower", "Higher")
# Bind columns together as dataset
results.table <- as.data.frame(cbind(HR, CI, P))
write.table(results.table,
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "coxfit_results_table_her-only.txt")





