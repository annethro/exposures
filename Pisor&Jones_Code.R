############# PISOR & JONES CODE #############

# Run with R version 4.0.2
# See GitHub for data dictionary/metadata: https://github.com/annethro/exposures
# To reduce participant identifiability, note that we converted participant ages to five-year bins after completing the modeling process. This does not affect model results (except for very small variations in parameter estimates), with the exception of the robustness check models in which reciprocity-based long-distance relationships were the outcome: in these models, you will notice that the credible intervals for extraversion do not include zero, despite what is reported in the main text. Note that the RIDs used in this data set also do not match the RIDs from other data sets in papers by AP, so these data cannot be merged.


### Load relevant libraries (install if needed using e.g., install.packages("mice") )

pckgs <- c("mice", "brms", "flextable", "officer", "data.table", "ggplot2")
lapply(pckgs, require, character.only = TRUE)


### Set to the working directory that contains the data file and load the data

setwd ("")

dat <- read.csv("Pisor&Jones_Exposures_Data.csv", header = T, stringsAsFactors = F) # Note missing data for one participant (RID 10) on extraversion measures. One option is to impute for them, as we do below; another option is to delete their data. As there is only one participant for whom data are missing, and only for a control variable, we opt to impute.


### "Center" (shift so that smallest value is zero) for P1-3
dat$D_Mean_Ctr <- dat$D_Mean - min(dat$D_Mean)
dat$E_Mean_Ctr <- dat$E_Mean - min(dat$E_Mean)
dat$Inter_Mean_Ctr <- dat$Inter_Mean - min(dat$Inter_Mean)
dat$D_Freq_Ctr <- dat$D_Freq - min(dat$D_Freq)
dat$E_Freq_Ctr <- dat$E_Freq - min(dat$E_Freq)


### Round participants who named three or more long-distance connections (n=27) down to 3 for comparable size across levels of LDRs named. Same for reciprocity-based long-distance relationships: n=21 named two or more, so bin for comparability.

dat$LDR <- ifelse (dat$LDR > 3, 3, dat$LDR)
dat$RLDR <- ifelse (dat$RLDR > 2, 2, dat$RLDR)


### Calculate participant age in 2017 (at the time of interview). Recall that YrBorn is binned into five-year categories; see note at the top of this file

dat$Age <- 2017 - dat$YrBorn - min (2017 - dat$YrBorn)


### Impute extraversion values for one participant missing data

## Tell R which variables (columns) to use in imputation for which variables (rows)

predMat <- matrix (rep (0, ncol (dat)^2), ncol = ncol (dat), nrow = ncol (dat))
fromthis <- ifelse (colnames (dat) %in% c("RID", "YrBorn", "Ex_Stranger", "Ex_Convo"), 0, 1) # In imputation using predictive mean matching, what happens is that the mice function compares the data of the person with missing values (RID 10) with other participants' data and fills in RID 10's missing values based on the match of their data with other people's data. Here, we are telling R not to impute extraversion from other participants' RIDs (that wouldn't make sense), from YrBorn because we also have Age in the data set and the two are perfectly collinear, and from Ex_Stranger and Ex_Convo because RID 10 doesn't HAVE data for those two variables.
predMat[ colnames (dat) %in% "Ex_Stranger", ] <- fromthis
predMat[ colnames (dat) %in% "Ex_Convo", ] <- fromthis

## Run imputation (and set seed for replicability; if you don't set a seed, R picks a random one each time and the imputed values move around)

dat1 <- complete (mice (dat, method = "pmm", predictorMatrix = predMat, seed = 1591413214)) # Seed chosen randomly with Sys.time() on June 5-20. PMM means predictive mean matching. Complete returns the data set (saved as dat1) with missing values populated with predicted values.


### Calculate summary measure of extraversion and "center"

dat1$Extra <- dat1$Ex_Convo + dat1$Ex_Stranger - min (dat1$Ex_Convo + dat1$Ex_Stranger)


### Remove variables not to be used in coming analyses

dat2 <- dat1 [ , ! (colnames(dat1) %in% c("YrBorn", "D_Mean", "E_Mean", "Inter_Mean", "Ex_Stranger", "Ex_Convo"))]

### Give LDR and RLDR appropriate class (integer) for running cratio models (which expect the outcome to be an integer, not numeric)
dat2$LDR <- ordered(dat2$LDR); dat2$RLDR <- ordered(dat2$RLDR)

### Set weakly informative priors for all models
priors <- c(set_prior("normal(0,10)", class = "b"),
            set_prior("cauchy(0,2)", class = "sd"))

############# HYPOTHESIS 1 #############

# Note that if you are using R 4.0 or later, even if disc (that is, discrimination; see Item Response Theory literature) is not specified in your model (e.g., p1d_ei), brms will ask Stan (via rstan) to estimate it. rstan will return an error that estimated sample sizes were too low and chains have not mixed. However, when you look at plots (e.g., plot(p1d_fio)) and summaries, you'll see that everything has mixed -- just that there is an empty plot for disc and that it has not been estimated (it's just a vector of 1s with a standard error of 0). At the time of writing, rstan had not remedied this issue. However, it does not negatively impact model fitting -- just provides an unfortunate source of confusion for the user.

# Reminder: your results will not be identical to those reported in the manuscript, as age has been binned into five year categories.


### Set a seed such that results are reproducible (NB: if you don't do this, estimates move around, but only by a little)
semilla <- 1593839097


####### PREDICTION 1 #######

### Longer average drought interval or longer average excess precipitation interval -> more long-distance relationships

##### DROUGHT ######

### Drought (d), thresholds as flexible (f), include imputed data (i), include outliers for P1 (o).
# We wished to assess equal variance between communities at the same time as checking the influence of outliers here, but model fit issues precluded this. Accordingly, we started with just looking at outliers. (Section 2.4.4)

p1d_fio <- brm(formula = bf(LDR ~ D_Mean_Ctr + 
      Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
      (1|Community)) + # Main model
        lf(disc ~ D_Mean_Ctr), # This estimation is first step in checking for influential outliers: it estimates the effect of each observation for D_Mean_Ctr on the residual standard deviation of the outcome.
    family = cratio (link = "logit", link_disc = "log", threshold = "flexible"), # Will compare to equidistance assumptions to see if those are warranted; they require estimating fewer parameters and so are preferable by principles of parsimony. The choice of a logit link for the main model and a log link for the residual standard deviation model are standard and well-supported for the models we're running here; see https://journals.sagepub.com/doi/10.1177/2515245918823199 for details.
    prior = priors, control = list(adapt_delta = 0.99, max_treedepth = 12), # You will see that we increased these throughout to help with model fitting. Adapt_delta decreases the step size -- the probability that the MCMC model fitting process accepts the next proposed posterior distribution (see https://link.springer.com/article/10.3758/s13428-016-0746-9 for an explanation). This increases the robustness of estimates of the posterior distribution. You don't need to normally change adapt_delta unless brm throws an error about lack of convergence; our models threw such errors. Max_treedepth just sets a limit to the number of trees (the values visited by a random walk in the fitting process before the sampler turns back; see https://arxiv.org/pdf/1111.4246.pdf) assessed by the modeling process, which keeps the model from running for an inordinately long time.
    data = dat2, chains = 4, seed = semilla) 

## Check whether outliers influence model fit.

plot(predict(p1d_fio, dpar = "disc")) # RIDs 67 and 75 are outliers for mean drought length, as demonstrated both by exploratory data analysis and plotting predicted values. Refit model with these removed.


### Drought, thresholds as flexible, include imputed data, NO outliers included, assess whether unequal variance by community (u) (see Supporting Information 3.3)

p1d_fiu <- brm(formula = bf(LDR ~ D_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) + # Main model
                 lf(disc ~ Community), # Assess whether unequal variance by community: how much does being in one community vs the other affect the residual standard deviation of the outcome?
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2[!(dat2$RID %in% c(67,75)),], chains = 4, seed = semilla) # Note that the two outliers are excluded here and in all other p1d (that is, prediction 1: drought) models.

# Use summary() to look at the model estimates. See how the CI for the Moseten community include 0. Accordingly, no evidence that unequal variance by community. (See SI Section 3.3)


### Drought, thresholds as flexible, include imputed data, NO outliers included

p1d_fi <- brm(formula = bf(LDR ~ D_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), # Main model
              family = cratio (link = "logit", threshold = "flexible"),
              prior = priors, control = list(adapt_delta = 0.99, max_treedepth = 15),
              data = dat2[!(dat2$RID %in% c(67,75)),], chains = 4, seed = semilla)


### Drought, thresholds as equidistant, include imputed data, NO outliers included
# This model requires estimating fewer parameters than the one immediately above. Accordingly, by the principles of parsimony, if a comparison of model fits suggests a better fit for this one, we will equidistant thresholds to report our results.

p1d_ei <- brm(formula = bf(LDR ~ D_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), # Main model
              family = cratio (link = "logit", threshold = "equidistant"),
              prior = priors, control = list(adapt_delta = 0.999, max_treedepth = 15),
              data = dat2[!(dat2$RID %in% c(67,75)),], chains = 4, seed = semilla)


### Both p1d_ei and p1d_fi show no issues with fit, although ei has more outliers in the random effect for community (something we assessed using a pairs() plot). Compare the two models to see if the assumption of equidistance is warranted.

loo(p1d_ei,p1d_fi)

# Leave-one-out cross-validation suggests that ei has a slightly better fit -- it has a slightly lower information criterion (LOOIC value; SI Section 3.3). We will evaluate model p1d_ei in the main text.


### Re-run ei with participant for whom we imputed extraversion data (RID 10) excluded as a robustness check (Section 2.4.5).

p1d_e <- brm(formula = bf(LDR ~ D_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), # Main model
              family = cratio (link = "logit", threshold = "equidistant"),
              prior = priors, control = list(adapt_delta = 0.99, max_treedepth = 15),
              data = dat2[!(dat2$RID %in% c(67,75,10)),], chains = 4, seed = semilla)

# Robust to their inclusion or exclusion (note that even though intercept estimates move, delta (the spacing of intervals) stays roughly the same).


### Re-run ei with robustness check of outcome: reciprocal long-distance relationships only (p1r). Use same model structure (ei) for comparability and check fit. (Section 2.4.3; to see model results, see Figure S1)

p1rd_ei <- brm(formula = bf(RLDR ~ D_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), # Main model
              family = cratio (link = "logit", threshold = "equidistant"),
              prior = priors, control = list(adapt_delta = 0.999, max_treedepth = 15),
              data = dat2[!(dat2$RID %in% c(67,75)),], chains = 4, seed = semilla)


##### EXCESS PRECIPITATON #####

# Same process as for drought for consistency.


### Excess precip, thresholds as flexible, include imputed data, include outliers for P1.
# We wished to assess equal variance between communities at the same time as checking the influence of outliers here, but model fit issues precluded this. Accordingly, we started with just looking at outliers. For more details on what everything means, see p1d_fio, above.

p1e_fio <- brm(formula = bf(LDR ~ E_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) + # Main model
                 lf(disc ~
                      E_Mean_Ctr), # This estimation is first step in checking for influential outliers: it estimates the effects of each variable on the residual standard deviation of the outcome.
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"), # Will compare to equidistance assumptions to see if those are warranted; they require estimating fewer parameters.
               prior = priors, control = list(adapt_delta = 0.99, max_treedepth = 12),
               data = dat2, chains = 4, seed = semilla) 

plot(predict(p1e_fio, dpar = "disc")) # No obvious outliers cross-checking this with exploratory data analysis (that is, where we just do plot(dat2$E_Mean_Ctr)).


### Excess precip, thresholds as flexible, include imputed data, assess whether unequal variance by community, include all data points (candidate outliers still in, as didn't influence fit, but to avoid confusion, we stop including the "o" in the model names)

p1e_fiu <- brm(formula = bf(LDR ~ E_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) + # Main model
                 lf(disc ~ Community), # Assess whether unequal variance by community
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2, chains = 4, seed = semilla) 

# No evidence of unequal variance by community.


### Excess precip, thresholds as flexible, include imputed data

p1e_fi <- brm(formula = bf(LDR ~ E_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)),
              family = cratio (link = "logit", threshold = "flexible"),
              prior = priors, control = list(adapt_delta = 0.99),
              data = dat2, chains = 4, seed = semilla)

#Made thresholds smaller. Compare to equidistant fit to see if equidistance assumptions are warranted.


### Excess precipitation, thresholds as equidistant, include imputed data

p1e_ei <- brm(formula = bf(LDR ~ E_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), # Main model
              family = cratio (link = "logit", threshold = "equidistant"),
              prior = priors, control = list(adapt_delta = 0.999, max_treedepth = 12),
              data = dat2, chains = 4, seed = semilla)

loo(p1e_ei,p1e_fi)

# Leave-one-out cross-validation suggests that ei has a slightly better fit -- it has a slightly lower information criterion (LOOIC value). We will evaluate model p1e_ei in the main text.


### Re-run ei with participant for whom we imputed excluded as a robustness check.

p1e_e <- brm(formula = bf(LDR ~ E_Mean_Ctr + 
                            Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                            (1|Community)), # Main model
             family = cratio (link = "logit", threshold = "equidistant"),
             prior = priors, control = list(adapt_delta = 0.999),
             data = dat2[!(dat2$RID %in% c(10)),], chains = 4, seed = semilla)

# Robust to everything EXCEPT sex falls from significance with this individual excluded (Section 3.5). However, that's not extraversion, which stays the same in its ("significant") effect. There was also nothing special about this person that led them to not complete the extraversion component: AP just forgot to do it. Leave this person in.


### Re-run ei with robustness check of outcome: reciprocal long-distance relationships only (p1r). Use same model structure (ei) for comparability and check fit.

p1re_ei <- brm(formula = bf(RLDR ~ E_Mean_Ctr + 
                            Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                            (1|Community)), # Main model
             family = cratio (link = "logit", threshold = "equidistant"),
             prior = priors, control = list(adapt_delta = 0.99),
             data = dat2, chains = 4, seed = semilla)


####### ALTERNATIVE PREDICTION 1 #######

### Longer average drought interval or longer average excess precipitation interval -> presence of same-community relationships (Section 2.4.3). Note that this is a binomial logistic regression, given the bernoulli (0,1) outcome.

# For comparability to other P1 models (above: p1d_ei, p1e_ei), use ei model formulation and just check fit.


##### DROUGHT ######

### Alternative outcome (p1a), drought (d), include imputed data (i) (see Figure S2 for model results)

p1ad_ei <- brm(formula = bf(SCR ~ D_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), # Main model
              family = bernoulli,
              prior = priors, control = list(adapt_delta = 0.999),
              data = dat2[!(dat2$RID %in% c(67,75)),], chains = 4, seed = semilla)

# A couple of potential outliers for community random effect (according to pairs plot) don't appear to be a problem when we plot the predicted values (plot(predict(p1ad_ei,dpar = "sd"))). Just a higher adapt_delta needed to avoid divergent transitions.


### Alternative outcome (p1a), drought (d), remove imputed data as a robustness check

p1ad_e <- brm(formula = bf(SCR ~ D_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), # Main model
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2[!(dat2$RID %in% c(67,75,10)),], chains = 4, seed = semilla)

# No change relative to ei.


### Alternative outcome: non-consanguineal kin only (p1an), drought (d), include imputed data (i) (Section 2.4.3; see Figure S3 for model results)

p1and_ei <- brm(formula = bf(KSCR ~ D_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), # Main model
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2[!(dat2$RID %in% c(67,75)),], chains = 4, seed = semilla)


##### EXCESS PRECIPITATION ######

### Alternative outcome (p1a), excess precipitation (e), include imputed data (i)

p1ae_ei <- brm(formula = bf(SCR ~ E_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), # Main model
              family = bernoulli,
              prior = priors, control = list(adapt_delta = 0.9999),
              data = dat2, chains = 4, seed = semilla)

# Again, pairs plot shows some candidate outliers for the random effect, but examining predicted values for the random effect suggest that none are influencing model fit.


### Alternative outcome (p1a), excess precipitation (e), remove imputed data as a robustness check

p1ae_e <- brm(formula = bf(SCR ~ E_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), # Main model
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2[!(dat2$RID %in% c(10)),], chains = 4, seed = semilla)

# Consistent


### Alternative outcome: non-consanguineal kin only (p1an), excess precipitation, include imputed data (i)

p1ane_ei <- brm(formula = bf(KSCR ~ E_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), # Main model
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2, chains = 4, seed = semilla)


############# HYPOTHESIS 2 #############

####### PREDICTION 2 #######

### Shorter average interval free of droughts or excess precipitation -> more long-distance relationships

# Reminder: there is only one set of models for this prediction, instead of separate models for drought and excess precipitation


### Thresholds as flexible (f), include imputed data (i), include outliers for P2 (o).

p2_fio <- brm(formula = bf(LDR ~ Inter_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) +
                 lf(disc ~ Inter_Mean_Ctr),
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.999, max_treedepth = 15),
               data = dat2, chains = 4, seed = semilla)

plot(predict(p2_fio, dpar = "disc")) # No likely influential points here. Moving on without removing participants.


### Thresholds as flexible, include imputed data, assess whether unequal variance by community, include all data points (candidate outliers still in, as didn't influence fit, but to avoid confusion, we stop including the "o" in the model names)

p2_fiu <- brm(formula = bf(LDR ~ Inter_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) +
                 lf(disc ~ Community), # Assess whether unequal variance by community
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2, chains = 4, seed = semilla) 

# No evidence for unequal variance by community.


### Thresholds as flexible, include imputed data

p2_fi <- brm(formula = bf(LDR ~ Inter_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)),
              family = cratio (link = "logit", threshold = "flexible"),
              prior = priors, control = list(adapt_delta = 0.99),
              data = dat2, chains = 4, seed = semilla)

# Compare to equidistant fit to see if equidistance assumptions are warranted.


### Excess precipitation, thresholds as equidistant, include imputed data

p2_ei <- brm(formula = bf(LDR ~ Inter_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), 
              family = cratio (link = "logit", threshold = "equidistant"),
              prior = priors, control = list(adapt_delta = 0.9999, max_treedepth = 12),
              data = dat2, chains = 4, seed = semilla)

loo(p2_ei,p2_fi)

# Leave-one-out cross-validation suggests that ei has a slightly better fit -- it has a slightly lower information criterion (LOOIC value). We will evaluate model p2_ei in the main text.


### Re-run ei with participant for whom we imputed excluded as a robustness check.

p2_e <- brm(formula = bf(LDR ~ Inter_Mean_Ctr + 
                            Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                            (1|Community)), 
             family = cratio (link = "logit", threshold = "equidistant"),
             prior = priors, control = list(adapt_delta = 0.9999, max_treedepth = 12),
             data = dat2[!(dat2$RID %in% c(10)),], chains = 4, seed = semilla)

# Robust.


### Re-run ei with robustness check of outcome: reciprocal long-distance relationships only (p2r). Use same model structure (ei) for comparability and check fit.

p2r_ei <- brm(formula = bf(RLDR ~ Inter_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), 
               family = cratio (link = "logit", threshold = "equidistant"),
               prior = priors, control = list(adapt_delta = 0.9999, max_treedepth = 12),
               data = dat2, chains = 4, seed = semilla)


####### ALTERNATIVE PREDICTION 2 #######

### Shorter in-between interval -> presence of same-community relationships. Note that this is a binomial logistic regression, given the bernoulli (0,1) outcome.

# For comparability to other P2 models, use ei model formulation and just check fit.


### Alternative outcome (p2a), include imputed data (i)

p2a_ei <- brm(formula = bf(SCR ~ Inter_Mean_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), 
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2, chains = 4, seed = semilla)


### Alternative outcome (p2a), remove imputed data as a robustness check

p2a_e <- brm(formula = bf(SCR ~ Inter_Mean_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), 
              family = bernoulli,
              prior = priors, control = list(adapt_delta = 0.999),
              data = dat2[!(dat2$RID %in% c(10)),], chains = 4, seed = semilla)

# No change relative to ei.


### Alternative outcome: non-consanguineal kin only (p2an), include imputed data (i)

p2an_ei <- brm(formula = bf(KSCR ~ Inter_Mean_Ctr + 
                               Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                               (1|Community)), 
                family = bernoulli,
                prior = priors, control = list(adapt_delta = 0.999),
                data = dat2, chains = 4, seed = semilla)


############# HYPOTHESIS 3 #############

# The less rare a downside (drought or excess precipitation) -> the more long-distance friends individuals should have.

####### PREDICTION 3 #######

##### DROUGHT ######

### Drought (d), thresholds as flexible (f), include imputed data (i), include outliers for P3 (o).

p3d_fio <- brm(formula = bf(LDR ~ D_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) +
                 lf(disc ~ D_Freq_Ctr), 
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2, chains = 4, seed = semilla)

## Check whether outliers influence model fit.

plot(predict(p3d_fio, dpar = "disc")) # No clear evidence of influence, and exclusion of one candidate (RID 58) reduces effective sample size. Don't omit anyone.


### Drought, thresholds as flexible, include imputed data, assess whether unequal variance by community (u)

p3d_fiu <- brm(formula = bf(LDR ~ D_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) +
                 lf(disc ~ Community), 
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2, chains = 4, seed = semilla) 

# No evidence of unequal variance by community.


### Drought, thresholds as flexible, include imputed data

p3d_fi <- brm(formula = bf(LDR ~ D_Freq_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)),
              family = cratio (link = "logit", threshold = "flexible"),
              prior = priors, control = list(adapt_delta = 0.99),
              data = dat2, chains = 4, seed = semilla)


### Drought, thresholds as equidistant, include imputed data

p3d_ei <- brm(formula = bf(LDR ~ D_Freq_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)),
              family = cratio (link = "logit", threshold = "equidistant"),
              prior = priors, control = list(adapt_delta = 0.999),
              data = dat2, chains = 4, seed = semilla)


### Both p3d_ei and p3d_fi show no issues with fit. Compare the two models to see if the assumption of equidistance is warranted.

loo(p3d_fi,p3d_ei)

# Leave-one-out cross-validation suggests that ei has a slightly better fit -- it has a slightly lower information criterion (LOOIC value). We will evaluate model p1d_ei in the main text.


### Re-run ei with participant for whom we imputed excluded as a robustness check.

p3d_e <- brm(formula = bf(LDR ~ D_Freq_Ctr + 
                            Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                            (1|Community)),
             family = cratio (link = "logit", threshold = "equidistant"),
             prior = priors, control = list(adapt_delta = 0.999),
             data = dat2[dat2$RID!=10,], chains = 4, seed = semilla)

#Robust to their inclusion or exclusion.


### Re-run ei with robustness check of outcome: reciprocal long-distance relationships only (p3r). Use same model structure (ei) for comparability and check fit.

p3rd_ei <- brm(formula = bf(RLDR ~ D_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), 
               family = cratio (link = "logit", threshold = "equidistant"),
               prior = priors, control = list(adapt_delta = 0.999, max_treedepth = 12),
               data = dat2, chains = 4, seed = semilla)


##### EXCESS PRECIPITATION ######

### Excess precip (e), thresholds as flexible (f), include imputed data (i), include outliers for P3 (o).

p3e_fio <- brm(formula = bf(LDR ~ E_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) +
                 lf(disc ~ E_Freq_Ctr), 
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2, chains = 4, seed = semilla)

## Check whether outliers influence model fit.

plot(predict(p3e_fi, dpar = "disc")) # Maybe RID 67, however, re-running without them doesn't much improve fit and just makes other points look like outliers. Keep in.


### Excess, thresholds as flexible, include imputed data, assess whether unequal variance by community (u)

p3e_fiu <- brm(formula = bf(LDR ~ E_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)) +
                 lf(disc ~ Community), 
               family = cratio (link = "logit", link_disc = "log", threshold = "flexible"),
               prior = priors, control = list(adapt_delta = 0.99),
               data = dat2, chains = 4, seed = semilla)

# No evidence of unequal variance by community.


### Drought, thresholds as flexible, include imputed data

p3e_fi <- brm(formula = bf(LDR ~ E_Freq_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)),
              family = cratio (link = "logit", threshold = "flexible"),
              prior = priors, control = list(adapt_delta = 0.999),
              data = dat2, chains = 4, seed = semilla)


### Drought, thresholds as equidistant, include imputed data

p3e_ei <- brm(formula = bf(LDR ~ E_Freq_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)),
              family = cratio (link = "logit", threshold = "equidistant"),
              prior = priors, control = list(adapt_delta = 0.999),
              data = dat2, chains = 4, seed = semilla)


### Both p3d_fi and p3d_ei show no issues with fit. Compare the two models to see if the assumption of equidistance is warranted.

loo(p3e_ei,p3e_fi)

# Leave-one-out cross-validation suggests that ei has a slightly better fit -- it has a slightly lower information criterion (LOOIC value). We will evaluate model p1d_ei in the main text.


### Re-run ei with participant for whom we imputed excluded as a robustness check.

p3e_e <- brm(formula = bf(LDR ~ E_Freq_Ctr + 
                            Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                            (1|Community)),
             family = cratio (link = "logit", threshold = "equidistant"),
             prior = priors, control = list(adapt_delta = 0.999),
             data = dat2[dat2$RID!=10,], chains = 4, seed = semilla)

#Robust to their inclusion or exclusion.


### Re-run ei with robustness check of outcome: reciprocal long-distance relationships only (p3r). Use same model structure (ei) for comparability and check fit.

p3re_ei <- brm(formula = bf(RLDR ~ E_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), 
               family = cratio (link = "logit", threshold = "equidistant"),
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2, chains = 4, seed = semilla)


####### ALTERNATIVE PREDICTION 3 #######

### The less rare a downside (drought or excess precipitation) -> the more likely a participant will be to have a same-community relationship.

# For comparability to other P3 models, use ei model formulation and just check fit.


##### DROUGHT ######

### Alternative outcome (p3a), drought (d), include imputed data (i)

p3ad_ei <- brm(formula = bf(SCR ~ D_Freq_Ctr +
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), 
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2, chains = 4, seed = semilla)


### Alternative outcome (p13, drought (d), remove imputed data as a robustness check

p3ad_e <- brm(formula = bf(SCR ~ D_Freq_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), 
              family = bernoulli,
              prior = priors, control = list(adapt_delta = 0.999),
              data = dat2[!(dat2$RID %in% c(10)),], chains = 4, seed = semilla)

# No change relative to ei.


### Alternative outcome: non-consanguineal kin only (p1an), drought (d), include imputed data (i)

p3and_ei <- brm(formula = bf(KSCR ~ D_Freq_Ctr + 
                               Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                               (1|Community)), 
                family = bernoulli,
                prior = priors, control = list(adapt_delta = 0.99),
                data = dat2, chains = 4, seed = semilla)


##### EXCESS PRECIPITATION ######

### Alternative outcome (p3a), excess precipitation (e), include imputed data (i)

p3ae_ei <- brm(formula = bf(SCR ~ E_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), 
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2, chains = 4, seed = semilla)

# Again, pairs plot shows some candidate outliers for the random effect, but examining predicted values for the random effect suggest that none are influencing model fit.


### Alternative outcome (p3a), excess precipitation (e), remove imputed data as a robustness check

p3ae_e <- brm(formula = bf(SCR ~ E_Freq_Ctr + 
                             Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                             (1|Community)), 
              family = bernoulli,
              prior = priors, control = list(adapt_delta = 0.99),
              data = dat2[!(dat2$RID %in% c(10)),], chains = 4, seed = semilla)

# Consistent


### Alternative outcome: non-consanguineal kin only (p3an), excess precipitation, include imputed data (i)

p3ane_ei <- brm(formula = bf(KSCR ~ E_Freq_Ctr + 
                              Sex + Age + Smartphone + Vehicle + Extra + Mo_Travel +
                              (1|Community)), 
               family = bernoulli,
               prior = priors, control = list(adapt_delta = 0.999),
               data = dat2, chains = 4, seed = semilla)


################## SUMMARY STATISTICS (TABLE S1) ###################

##### Ordinal table, levels 0-3 #####

zerothru <- list("LDR" = summary(as.factor(dat$LDR)), "RLDR" = summary(as.factor(dat$RLDR)), "SCR" = summary(as.factor(dat$SCR)), "KSCR" = summary(as.factor(dat$KSCR)), "Smartphone" = summary(as.factor(dat$Smartphone)), "Vehicle" = summary(as.factor(dat$Vehicle)))

zerothru_tab <- data.frame(t(data.frame(lapply(zerothru, "length<-", max(lengths(zerothru)))))) # Makes a table with counts for each level of each variable
zerothru_tab[is.na(zerothru_tab)] <- ""
setDT(zerothru_tab, keep.rownames = TRUE)[] # Convert the rownames from the data frame into an actual column, necessary for the flextable command below.

colnames(zerothru_tab)[colnames(zerothru_tab) %in% "rn"] <- "Variable"
colnames(zerothru_tab)[colnames(zerothru_tab) %in% "X0"] <- "Zero"
colnames(zerothru_tab)[colnames(zerothru_tab) %in% "X1"] <- "One"
colnames(zerothru_tab)[colnames(zerothru_tab) %in% "X2"] <- "Two"
colnames(zerothru_tab)[colnames(zerothru_tab) %in% "X3"] <- "Three"

zerothru_tab$Variable[zerothru_tab$Variable=="LDR"] <- "Long-Distance"
zerothru_tab$Variable[zerothru_tab$Variable=="RLDR"] <- "Reciprocal Long-Distance"
zerothru_tab$Variable[zerothru_tab$Variable=="SCR"] <- "Same-Community"
zerothru_tab$Variable[zerothru_tab$Variable=="KSCR"] <- "Non-Con. Kin Same-Community"

zt <- flextable(zerothru_tab) # flextable creates good-looking tables that can be output in R Markdown, Latex, and Office.
zt <- align_text_col(zt, align = "center", header = TRUE)
zt <- align_nottext_col(zt, align = "center", header = TRUE)
zt <- width(zt, j = "Variable", width = 2.5)


##### Ordinal table, 1-4 #####

onethru <- list("Ex_Stranger" = summary(factor(dat$Ex_Stranger[!is.na(dat$Ex_Stranger)])), "Ex_Convo" = summary(as.factor(dat$Ex_Convo[!is.na(dat$Ex_Convo)])))

onethru_tab <- data.frame(t(data.frame(lapply(onethru, "length<-", max(lengths(onethru))))))
onethru_tab[is.na(onethru_tab)] <- ""
setDT(onethru_tab, keep.rownames = TRUE)[]

colnames(onethru_tab)[colnames(onethru_tab) %in% "rn"] <- "Variable"
colnames(onethru_tab)[colnames(onethru_tab) %in% "X1"] <- "One"
colnames(onethru_tab)[colnames(onethru_tab) %in% "X2"] <- "Two"
colnames(onethru_tab)[colnames(onethru_tab) %in% "X3"] <- "Three"
colnames(onethru_tab)[colnames(onethru_tab) %in% "V4"] <- "Four"

onethru_tab$Variable[onethru_tab$Variable=="Ex_Stranger"] <- "Extraversion: Stranger"
onethru_tab$Variable[onethru_tab$Variable=="Ex_Convo"] <- "Extraversion: Conversation"

ot <- flextable(onethru_tab)
ot <- align_text_col(ot, align = "center", header = TRUE)
ot <- align_nottext_col(ot, align = "center", header = TRUE)
ot <- width(ot, j = "Variable", width = 2)


##### Nominal table #####
  
nominal <-  data.frame(Variable= c("Community", "Sex"), rbind(as.numeric(summary(as.factor(dat$Community))), as.numeric(summary(as.factor(dat$Sex)))))
colnames(nominal)[colnames(nominal) %in% c("X1", "X2")] <- c("One", "Two")

nt <- flextable(nominal)
nt <- align_text_col(nt, align = "center", header = TRUE)
nt <- align_nottext_col(nt, align = "center", header = TRUE)
nt <- width(nt, j = "Variable", width = 1)


##### Continuous table #####

cont_mn <- c("D_Freq" = mean(dat$D_Freq), "E_Freq" = mean(dat$E_Freq), "D_Mean" = mean(dat$D_Mean), "E_Mean" = mean(dat$E_Mean), "Inter_Mean" = mean(dat$Inter_Mean), "Yr_Born" = mean(dat$YrBorn), "Mo_Travel" = mean(dat$Mo_Travel))

cont_sd <- c("D_Freq" = sd(dat$D_Freq), "E_Freq" = sd(dat$E_Freq), "D_Mean" = sd(dat$D_Mean), "E_Mean" = sd(dat$E_Mean), "Inter_Mean" = sd(dat$Inter_Mean), "Yr_Born" = sd(dat$YrBorn), "Mo_Travel" = sd(dat$Mo_Travel))

cont_tab <- data.frame(Mean = cont_mn, SD = cont_sd) # Rather than counts by level, for continuous variables, create a table of means and standard deviations for each.
setDT(cont_tab, keep.rownames = TRUE)[]
colnames(cont_tab)[colnames(cont_tab) %in% "rn"] <- "Variable"
cont_tab[ ,c("Mean", "SD")] <- lapply(cont_tab[ ,c("Mean", "SD")], function(x) {round(x,2)}) # Round the values on the table to make it less unwiedly.

cont_tab$Variable[cont_tab$Variable %in% "D_Freq"] <- "Percent Dry Months"
cont_tab$Variable[cont_tab$Variable %in% "E_Freq"] <- "Percent Wet Months"
cont_tab$Variable[cont_tab$Variable %in% "D_Mean"] <- "Mean Len. Drought"
cont_tab$Variable[cont_tab$Variable %in% "E_Mean"] <- "Mean Len. Excess P."
cont_tab$Variable[cont_tab$Variable %in% "Inter_Mean"] <- "Mean Len. No D. or E.P."
cont_tab$Variable[cont_tab$Variable %in% "Yr_Born"] <- "Birth Year"
cont_tab$Variable[cont_tab$Variable %in% "Mo_Travel"] <- "Depts. & Countries Visited"

ct <- flextable(cont_tab)
ct <- align_text_col(ct, align = "center", header = TRUE)
ct <- align_nottext_col(ct, align = "center", header = TRUE)
ct <- width(ct, j = "Variable", width = 2.5)

##### Export all tables #####

# Exporting all tables to a single Word document using the OfficeR package. Note that OfficeR can export to other Office platforms; to export to Latex or R Markdown, see documentation for the flextable package.

word_export <- read_docx()
body_add_flextable(word_export, zt)
body_add_par(word_export, value = "")
body_add_flextable(word_export, ot)
body_add_par(word_export, value = "")
body_add_flextable(word_export, nt)
body_add_par(word_export, value = "")
body_add_flextable(word_export, ct)
print(word_export, 'summary_stats.docx')


################## Figure illustrating the three predictions (Figure 1) ###################

### Pred 1: Long duration vs short duration

pdf("Figure 1a.pdf", 4, 4)

par(mar = c(1.5, 2.5, 0.6, 0.5), mgp = c(0.5, 0, 0))

plot(0, type = "n", xlim = c(-1,2), ylim = c(0,2), xaxt = "n", yaxt = "n", xlab = "Time", ylab = "Below-Mean Resources", main = "", cex.lab = 1.5)

text(-1.5, 2.08, labels = "A", xpd = NA, cex = 1.5, font = 2, family = "serif")

dens_a1 <- density(runif(10000000, 0, 1.25) - 1)
x_a1 <- dens_a1$x
y_a1 <- dens_a1$y

dens_a2 <- density(runif(10000000, 0, 0.5) + 1)
x1_a2 <- dens_a2$x
y1_a2 <- dens_a2$y

lines(x_a1[y_a1 <= 0.5], y_a1[y_a1 <= 0.5], lwd=5)
lines(x1_a2[y1_a2 <= 0.5], y1_a2[y1_a2 <= 0.5], lwd=5)

text(-0.4, 0.5, pos = 3, "Longer", cex = 1.5, family = "serif")
text(1.3, 0.5, pos = 3, "Shorter", cex = 1.5, family = "serif")

dev.off()

### Pred 2: More autocorrelated or less

pdf("Figure 1b.pdf", 4, 4)

par(mar = c(1.5, 2.5, 0.6, 0.5), mgp = c(0.5, 0, 0))

plot(0, type = "n", xlim = c(-1, 100), ylim = c(0, 0.5), xaxt = "n", yaxt = "n", xlab = "Time", ylab = "Below-Mean Resources", main = "", cex.lab = 1.5)

text(-18, 0.52, labels = "B", xpd = NA, cex = 1.5, font = 2, family = "serif")

dens_b <- density(rnorm(100000, 5, 2))
x_b <- dens_b$x
y_b <- dens_b$y

lines(x_b[x_b < 9.5], y_b[x_b < 9.5],lwd=5)
lines(x_b[x_b > 0.75] + 10, y_b [x_b > 0.75],lwd=5)
lines(x_b + 25, y_b,lwd=5)

lines(x_b + 60, y_b,lwd=5)
lines(x_b + 90, y_b,lwd=5)

text(x = 18, y = 0.24, '{', srt = 270, cex = 6, family = 'serif')
text(x = 80, y = 0.24, '{', srt = 270, cex = 6, family = 'serif')

text(25, 0.32, pos = 3, "More \n Autocorrelated", cex = 1.5, family = 'serif', srt = 45)
text(80, 0.32, pos = 3, "Less \n Autocorrelated", cex = 1.5, family = 'serif', srt = 45)

dev.off()

### Pred 3: More patterned or more like a shock

pdf("Figure 1c.pdf", 4, 4)

par(mar = c(1.5, 2.5, 0.6, 0.5), mgp = c(0.5, 0, 0))

plot(0, type = "n", xlim = c(0,10), ylim = c(0,4), xaxt = "n", yaxt = "n", xlab = "Time", ylab = "Below-Mean Resources", main = "", cex.lab = 1.5)

text(-1.65, 4.15, labels = "C", xpd = NA, cex = 1.5, font = 2, family = "serif")

dens_c <- density(rnorm(100000, 3, 0.3))
x_c <- dens_c$x
y_c <- dens_c$y

lines(density(rbeta(1000000, 1, 2) + 0.5), lwd = 5)
lines(x_c[x_c > 2 & x_c < 4], y_c[x_c > 2 & x_c < 4], lwd=5)
lines(density(rbeta(1000000, 2, 1) + 8), lwd = 5)

text(x = 2, y = 2.25, '{', srt = 270, cex = 6, family = 'serif')
text(x = 8, y = 2.25, '{', srt = 270, cex = 6, family = 'serif')

text(2.5, 2.7, pos=3,"More\nPatterned",cex=1.5, family = 'serif', srt = 45)
text(8.5, 2.7, pos=3,"Less\nPatterned",cex=1.5, family = 'serif', srt = 45)

dev.off()


################## A FOREST PLOT INSPIRED BY ECON-STYLE TABLES FOR REPORTING DIFFERENT MODEL FITS (FIGURE 2) ###################

# Models are p1d_ei, p1e_ei, p2_ei, p3d_ei, p3e_ei

est_p1d <- data.frame(fixef(p1d_ei, probs = c(0.05, 0.95))) # This pulls the parameter estimates for the fixed effects out of the model fit, including 90% credible intervals.
setDT(est_p1d, keep.rownames = TRUE)[] # Convert row names into a column for use with the flextable package.

est_p1e <- data.frame(fixef(p1e_ei, probs = c(0.05, 0.95)))
setDT(est_p1e, keep.rownames = TRUE)[]

est_p2 <- data.frame(fixef(p2_ei, probs = c(0.05, 0.95)))
setDT(est_p2, keep.rownames = TRUE)[]

est_p3d <- data.frame(fixef(p3d_ei, probs = c(0.05, 0.95)))
setDT(est_p3d, keep.rownames = TRUE)[]

est_p3e <- data.frame(fixef(p3e_ei, probs = c(0.05, 0.95)))
setDT(est_p3e, keep.rownames = TRUE)[]

### Make a data frame for each model fit, then smash them into a single data frame using rbind.

df_p1d <- data.frame(Variable = est_p1d$rn, Group = "P1 Drought", LI = exp(est_p1d$Q5), Est = exp(est_p1d$Estimate), HI = exp(est_p1d$Q95))
df_p1e <- data.frame(Variable = est_p1e$rn, Group = "P1 Excess", LI =  exp(est_p1e$Q5), Est =  exp(est_p1e$Estimate), HI =  exp(est_p1e$Q95))
df_p2 <- data.frame(Variable = est_p2$rn, Group = "P2", LI =  exp(est_p2$Q5), Est =  exp(est_p2$Estimate), HI =  exp(est_p2$Q95))
df_p3d <- data.frame(Variable = est_p3d$rn, Group = "P3 Drought", LI =  exp(est_p3d$Q5), Est =  exp(est_p3d$Estimate), HI =  exp(est_p3d$Q95))
df_p3e <- data.frame(Variable = est_p3e$rn, Group = "P3 Excess", LI =  exp(est_p3e$Q5), Est =  exp(est_p3e$Estimate), HI =  exp(est_p3e$Q95))


df <- rbind(df_p1d, df_p1e, df_p2, df_p3d, df_p3e)
df1 <- df[!(df$Variable %in% grep("Intercept", df$Variable, value = T)), ] # Intercepts are not plotted, particularly because these are the parameter estimates for the thresholds between each level in cratio models -- and these are notoriously difficult to interpret (see https://journals.sagepub.com/doi/10.1177/2515245918823199).
df1$Type <- ifelse(df1$Variable %in% c("D_Mean_Ctr", "E_Mean_Ctr", "Inter_Mean_Ctr", "D_Freq_Ctr", "E_Freq_Ctr"), "Predictor", "Third Variable")

df1$Variable[df1$Variable %in% "D_Freq_Ctr"] <- "Percent Dry Months"
df1$Variable[df1$Variable %in% "E_Freq_Ctr"] <- "Percent Wet Months"
df1$Variable[df1$Variable %in% "D_Mean_Ctr"] <- "Mean Length Drought"
df1$Variable[df1$Variable %in% "E_Mean_Ctr"] <- "Mean Length Excess"
df1$Variable[df1$Variable %in% "Inter_Mean_Ctr"] <- "Mean Length No D. or E."
df1$Variable[df1$Variable %in% "Mo_Travel"] <- "Depts. & Countries Visited"
df1$Variable[df1$Variable=="Extra"] <- "Extraversion"
df1$Variable[df1$Variable=="Sexmale"] <- "Sex: Male"

df1$Variable <- factor(df1$Variable, levels = c("Percent Wet Months", "Percent Dry Months", "Mean Length No D. or E.", "Mean Length Excess", "Mean Length Drought", "Depts. & Countries Visited", "Extraversion", "Vehicle", "Smartphone", "Sex: Male", "Age")) #This re-orders the data frame such that variables appear in the order discussed in the manuscript.


### Plot it

p <- ggplot(df1,aes(x = Variable, y = Est, ymin = LI, ymax = HI)) + 
  geom_hline(aes(yintercept = 1), color = "gray50", linetype = "dashed") +
  geom_linerange(size = 1) + 
  geom_point(size = 2) +
  facet_grid(Type ~ Group, scales = "free_y") + 
  theme(strip.text = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), panel.grid = element_blank()) + 
  theme(strip.text.y = element_text(angle = 360))  + 
  coord_flip() + theme(panel.spacing = unit(1, "lines")) + 
  theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid")) + 
  theme(strip.background = element_rect(fill = "white", color = "black", linetype = "solid")) +
  labs(y = "Odds Ratios (ORs)")

ggsave("Figure 2.pdf", p, height=8, width=14)


################## FIGURE S1 (RECIPROCAL) ##################

# Models are p1rd_ei, p1re_ei, p2r_ei, p3rd_ei, p3re_ei

est_p1d <- data.frame(fixef(p1rd_ei, probs = c(0.05, 0.95)))
setDT(est_p1d, keep.rownames = TRUE)[]

est_p1e <- data.frame(fixef(p1re_ei, probs = c(0.05, 0.95)))
setDT(est_p1e, keep.rownames = TRUE)[]

est_p2 <- data.frame(fixef(p2r_ei, probs = c(0.05, 0.95)))
setDT(est_p2, keep.rownames = TRUE)[]

est_p3d <- data.frame(fixef(p3rd_ei, probs = c(0.05, 0.95)))
setDT(est_p3d, keep.rownames = TRUE)[]

est_p3e <- data.frame(fixef(p3re_ei, probs = c(0.05, 0.95)))
setDT(est_p3e, keep.rownames = TRUE)[]

df_p1d <- data.frame(Variable = est_p1d$rn, Group = "P1 Drought", LI = exp(est_p1d$Q5), Est = exp(est_p1d$Estimate), HI = exp(est_p1d$Q95))
df_p1e <- data.frame(Variable = est_p1e$rn, Group = "P1 Excess", LI =  exp(est_p1e$Q5), Est =  exp(est_p1e$Estimate), HI =  exp(est_p1e$Q95))
df_p2 <- data.frame(Variable = est_p2$rn, Group = "P2", LI =  exp(est_p2$Q5), Est =  exp(est_p2$Estimate), HI =  exp(est_p2$Q95))
df_p3d <- data.frame(Variable = est_p3d$rn, Group = "P3 Drought", LI =  exp(est_p3d$Q5), Est =  exp(est_p3d$Estimate), HI =  exp(est_p3d$Q95))
df_p3e <- data.frame(Variable = est_p3e$rn, Group = "P3 Excess", LI =  exp(est_p3e$Q5), Est =  exp(est_p3e$Estimate), HI =  exp(est_p3e$Q95))

df <- rbind(df_p1d, df_p1e, df_p2, df_p3d, df_p3e)
df1 <- df[!(df$Variable %in% grep("Intercept", df$Variable, value = T)), ]
df1$Type <- ifelse(df1$Variable %in% c("D_Mean_Ctr", "E_Mean_Ctr", "Inter_Mean_Ctr", "D_Freq_Ctr", "E_Freq_Ctr"), "Predictor", "Third Variable")

df1$Variable[df1$Variable %in% "D_Freq_Ctr"] <- "Percent Dry Months"
df1$Variable[df1$Variable %in% "E_Freq_Ctr"] <- "Percent Wet Months"
df1$Variable[df1$Variable %in% "D_Mean_Ctr"] <- "Mean Length Drought"
df1$Variable[df1$Variable %in% "E_Mean_Ctr"] <- "Mean Length Excess"
df1$Variable[df1$Variable %in% "Inter_Mean_Ctr"] <- "Mean Length No D. or E."
df1$Variable[df1$Variable %in% "Mo_Travel"] <- "Depts. & Countries Visited"
df1$Variable[df1$Variable=="Extra"] <- "Extraversion"
df1$Variable[df1$Variable=="Sexmale"] <- "Sex: Male"

df1$Variable <- factor(df1$Variable, levels = c("Percent Wet Months", "Percent Dry Months", "Mean Length No D. or E.", "Mean Length Excess", "Mean Length Drought", "Depts. & Countries Visited", "Extraversion", "Vehicle", "Smartphone", "Sex: Male", "Age"))

p <- ggplot(df1,aes(x = Variable, y = Est, ymin = LI, ymax = HI)) + 
  geom_hline(aes(yintercept = 1), color = "gray50", linetype = "dashed") +
  geom_linerange(size = 1) + 
  geom_point(size = 2) +
  facet_grid(Type ~ Group, scales = "free_y") + 
  theme(strip.text = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), panel.grid = element_blank()) + 
  theme(strip.text.y = element_text(angle = 360))  + 
  coord_flip() + theme(panel.spacing = unit(1, "lines")) + 
  theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid")) + 
  theme(strip.background = element_rect(fill = "white", color = "black", linetype = "solid")) +
  labs(y = "Odds Ratios (ORs)")

ggsave("Figure S1.pdf", p, height=8, width=14)


################## FIGURE S2 (SCR) ##################

# Models are p1ad_ei, p1ae_ei, p2a_ei, p3ad_ei, p3ae_ei

est_p1d <- data.frame(fixef(p1ad_ei, probs = c(0.05, 0.95)))
setDT(est_p1d, keep.rownames = TRUE)[]

est_p1e <- data.frame(fixef(p1ae_ei, probs = c(0.05, 0.95)))
setDT(est_p1e, keep.rownames = TRUE)[]

est_p2 <- data.frame(fixef(p2a_ei, probs = c(0.05, 0.95)))
setDT(est_p2, keep.rownames = TRUE)[]

est_p3d <- data.frame(fixef(p3ad_ei, probs = c(0.05, 0.95)))
setDT(est_p3d, keep.rownames = TRUE)[]

est_p3e <- data.frame(fixef(p3ae_ei, probs = c(0.05, 0.95)))
setDT(est_p3e, keep.rownames = TRUE)[]

df_p1d <- data.frame(Variable = est_p1d$rn, Group = "P1 Drought", LI = exp(est_p1d$Q5), Est = exp(est_p1d$Estimate), HI = exp(est_p1d$Q95))
df_p1e <- data.frame(Variable = est_p1e$rn, Group = "P1 Excess", LI =  exp(est_p1e$Q5), Est =  exp(est_p1e$Estimate), HI =  exp(est_p1e$Q95))
df_p2 <- data.frame(Variable = est_p2$rn, Group = "P2", LI =  exp(est_p2$Q5), Est =  exp(est_p2$Estimate), HI =  exp(est_p2$Q95))
df_p3d <- data.frame(Variable = est_p3d$rn, Group = "P3 Drought", LI =  exp(est_p3d$Q5), Est =  exp(est_p3d$Estimate), HI =  exp(est_p3d$Q95))
df_p3e <- data.frame(Variable = est_p3e$rn, Group = "P3 Excess", LI =  exp(est_p3e$Q5), Est =  exp(est_p3e$Estimate), HI =  exp(est_p3e$Q95))

df <- rbind(df_p1d, df_p1e, df_p2, df_p3d, df_p3e)
df1 <- df[!(df$Variable %in% grep("Intercept", df$Variable, value = T)), ]
df1$Type <- ifelse(df1$Variable %in% c("D_Mean_Ctr", "E_Mean_Ctr", "Inter_Mean_Ctr", "D_Freq_Ctr", "E_Freq_Ctr"), "Predictor", "Third Variable")

df1$Variable[df1$Variable %in% "D_Freq_Ctr"] <- "Percent Dry Months"
df1$Variable[df1$Variable %in% "E_Freq_Ctr"] <- "Percent Wet Months"
df1$Variable[df1$Variable %in% "D_Mean_Ctr"] <- "Mean Length Drought"
df1$Variable[df1$Variable %in% "E_Mean_Ctr"] <- "Mean Length Excess"
df1$Variable[df1$Variable %in% "Inter_Mean_Ctr"] <- "Mean Length No D. or E."
df1$Variable[df1$Variable %in% "Mo_Travel"] <- "Depts. & Countries Visited"
df1$Variable[df1$Variable=="Extra"] <- "Extraversion"
df1$Variable[df1$Variable=="Sexmale"] <- "Sex: Male"

df1$Variable <- factor(df1$Variable, levels = c("Percent Wet Months", "Percent Dry Months", "Mean Length No D. or E.", "Mean Length Excess", "Mean Length Drought", "Depts. & Countries Visited", "Extraversion", "Vehicle", "Smartphone", "Sex: Male", "Age"))

p <- ggplot(df1,aes(x = Variable, y = Est, ymin = LI, ymax = HI)) + 
  geom_hline(aes(yintercept = 1), color = "gray50", linetype = "dashed") +
  geom_linerange(size = 1) + 
  geom_point(size = 2) +
  facet_grid(Type ~ Group, scales = "free_y") + 
  theme(strip.text = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), panel.grid = element_blank()) + 
  theme(strip.text.y = element_text(angle = 360))  + 
  coord_flip() + theme(panel.spacing = unit(1, "lines")) + 
  theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid")) + 
  theme(strip.background = element_rect(fill = "white", color = "black", linetype = "solid")) +
  labs(y = "Odds Ratios (ORs)")

ggsave("Figure S2.pdf", p, height=8, width=14)


################## FIGURE S3 (KSCR) ##################

# Models are p1and_ei, p1ane_ei, p2an_ei, p3and_ei, p3ane_ei

est_p1d <- data.frame(fixef(p1and_ei, probs = c(0.05, 0.95)))
setDT(est_p1d, keep.rownames = TRUE)[]

est_p1e <- data.frame(fixef(p1ane_ei, probs = c(0.05, 0.95)))
setDT(est_p1e, keep.rownames = TRUE)[]

est_p2 <- data.frame(fixef(p2an_ei, probs = c(0.05, 0.95)))
setDT(est_p2, keep.rownames = TRUE)[]

est_p3d <- data.frame(fixef(p3and_ei, probs = c(0.05, 0.95)))
setDT(est_p3d, keep.rownames = TRUE)[]

est_p3e <- data.frame(fixef(p3ane_ei, probs = c(0.05, 0.95)))
setDT(est_p3e, keep.rownames = TRUE)[]

df_p1d <- data.frame(Variable = est_p1d$rn, Group = "P1 Drought", LI = exp(est_p1d$Q5), Est = exp(est_p1d$Estimate), HI = exp(est_p1d$Q95))
df_p1e <- data.frame(Variable = est_p1e$rn, Group = "P1 Excess", LI =  exp(est_p1e$Q5), Est =  exp(est_p1e$Estimate), HI =  exp(est_p1e$Q95))
df_p2 <- data.frame(Variable = est_p2$rn, Group = "P2", LI =  exp(est_p2$Q5), Est =  exp(est_p2$Estimate), HI =  exp(est_p2$Q95))
df_p3d <- data.frame(Variable = est_p3d$rn, Group = "P3 Drought", LI =  exp(est_p3d$Q5), Est =  exp(est_p3d$Estimate), HI =  exp(est_p3d$Q95))
df_p3e <- data.frame(Variable = est_p3e$rn, Group = "P3 Excess", LI =  exp(est_p3e$Q5), Est =  exp(est_p3e$Estimate), HI =  exp(est_p3e$Q95))

df <- rbind(df_p1d, df_p1e, df_p2, df_p3d, df_p3e)
df1 <- df[!(df$Variable %in% grep("Intercept", df$Variable, value = T)), ]
df1$Type <- ifelse(df1$Variable %in% c("D_Mean_Ctr", "E_Mean_Ctr", "Inter_Mean_Ctr", "D_Freq_Ctr", "E_Freq_Ctr"), "Predictor", "Third Variable")

df1$Variable[df1$Variable %in% "D_Freq_Ctr"] <- "Percent Dry Months"
df1$Variable[df1$Variable %in% "E_Freq_Ctr"] <- "Percent Wet Months"
df1$Variable[df1$Variable %in% "D_Mean_Ctr"] <- "Mean Length Drought"
df1$Variable[df1$Variable %in% "E_Mean_Ctr"] <- "Mean Length Excess"
df1$Variable[df1$Variable %in% "Inter_Mean_Ctr"] <- "Mean Length No D. or E."
df1$Variable[df1$Variable %in% "Mo_Travel"] <- "Depts. & Countries Visited"
df1$Variable[df1$Variable=="Extra"] <- "Extraversion"
df1$Variable[df1$Variable=="Sexmale"] <- "Sex: Male"

df1$Variable <- factor(df1$Variable, levels = c("Percent Wet Months", "Percent Dry Months", "Mean Length No D. or E.", "Mean Length Excess", "Mean Length Drought", "Depts. & Countries Visited", "Extraversion", "Vehicle", "Smartphone", "Sex: Male", "Age"))

p <- ggplot(df1,aes(x = Variable, y = Est, ymin = LI, ymax = HI)) + 
  geom_hline(aes(yintercept = 1), color = "gray50", linetype = "dashed") +
  geom_linerange(size = 1) + 
  geom_point(size = 2) +
  facet_grid(Type ~ Group, scales = "free_y") + 
  theme(strip.text = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), panel.grid = element_blank()) + 
  theme(strip.text.y = element_text(angle = 360))  + 
  coord_flip() + theme(panel.spacing = unit(1, "lines")) + 
  theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid")) + 
  theme(strip.background = element_rect(fill = "white", color = "black", linetype = "solid")) +
  labs(y = "Odds Ratios (ORs)")

ggsave("Figure S3.pdf", p, height=8, width=14)
