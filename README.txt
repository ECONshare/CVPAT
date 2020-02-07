###########################################################
## Input to CVPAT function
###########################################################

MV: 
Data of the manifest variables (indicators) used for the partial least squares (PLS) path model estimation. The variable names must be the same as specified in Model1 and Model2 (see below).

CVFOLDS: 
Sets the number of cross-validation folds.

Model1: 
The model specifications as required by the matrixpls package. That is, Model1 is a list of three elements. The first element is a matrix of the inner model, the second element is a matrix of the reflective measurement model and the last element is a matrix of the formative measurement model.

Model2: 
The structure is the same as for Model1, but the specifications apply to Model2.

hypothesis: 
Specifies the hypothesis to be tested. If hypothesis = "M1_better_out_of_sample_than_M2", then CVPAT will test whether M1 has a significantly better performance out-of-sample performance than M2. If Hypothesis = "M1!=M2", then CVPAT will test if the out-of-sample performance between M1 and M2 differs significantly from each other.

BootSamp: 
Sets the number of bootstrap samples.

boot.Di: 
If TRUE, bootstrapping is performed on the losses from a single cross-validation. If FALSE, bootstrapping is applied to the manifest variables, and consequently the PLS estimation and cross-validation are performed for each bootstrapping run (this option is much slower than boot.Di=TRUE)

seed: 
Gives the possibility to make reproducible results using the CVPAT function. 
If FALSE, the random number generator seed will change each time the CVPAT function is run. 
If the seed = x (where x is an integer; e.g. seed = 42), the random number generator uses a fixed seed x within the CVPAT function. If a seed = x is used, the seed is reset to the original value of the global environment after the CVPAT function is executed.



###########################################################
## Output from CVPAT function
###########################################################
The output from the CVPAT function is a list. The list has following elements:
- boot.p.values
- losses
- t.stat
- p.value
- conf.int
- conv.fail

----------
The content of each list element is:
----------

boot.p.values:
- p.value.perc.t: p-value calculated from the bootstrapped percentiles of the t-statistics.
- p.value.b.v.t: p-value calculated from the original t-statistic but replacing the variance of D_bar with the bootstrapped variance of D_bar
- p.value.perc.D: p-value calculated from the bootstrapped percentiles of D_bar.

losses$case_wise_loss:
- LossM1: The loss for each case for Model1
- LossM2: The loss for each case for Model2
- LossM1_sepLV: The loss for each case for each endogenous construct for Model1
- LossM2_sepLV: The loss for each case on each endogenous construct for Model2

losses$avg_losses:
- avg_losses_M1: The average loss overall for M1 and the average loss for each endogenous construct in M1
- avg_losses_M2: The average loss overall for M2 and the average loss for each endogenous construct in M2

t.stat: 
- t.stat: The non-bootstrapped t-statistic
- t.stat.b.v: The t-statistic using the bootstrapped variance of D_bar

p.val: 
The non-bootstrapped p-value

conf.int: 
The non-bootstrapped confidence interval

conv.fail: 
The proportion of bootstrap runs where the PLS-SEM algorithm had a convergence error in one or more of the cross-validation runs. Only relevant if boot.Di = FALSE.
