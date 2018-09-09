# !diagnostics off
start_time <- Sys.time()


library(survAUC) ## For Schemper and Hendersons R2
library(rms) ## to fit Cox model with cph 



library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(survival)
library(parallel)
library(caret)
library(rsample)
library(VIM)
library(rstan)
library(rstanarm)
theme_set(theme_bw())
devtools::document()
library(survivalROC)

library(survbayes2)
options(mc.cores = parallel::detectCores() )
rstan_options(auto_write = TRUE)

data(lusc)

#convert variable names to lower case
colnames(lusc) <- tolower(colnames(lusc))
#convert all string NA to NA format
lusc <- lusc %>%
  mutate_all(funs(convert_blank_to_na))

lusc %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE,
            sortVars = TRUE, sortCombs = TRUE, plot = FALSE, only.miss = TRUE)

# Filter only smoking population
# lusc <- lusc %>%
#   dplyr::filter(tobacco_smoking_history_indicator != "1" )

lusc <- lusc %>% 
  dplyr::select( patient_id, os_status, os_months,  ajcc_nodes_pathologic_pn, ajcc_tumor_pathologic_pt, age,  icd_10, sex , ajcc_pathologic_tumor_stage, ajcc_metastasis_pathologic_pm, tobacco_smoking_history_indicator ) %>%
  dplyr::mutate(age = as.numeric(age),
                os_months = as.numeric(os_months),
                os_status = as.logical(os_status == "DECEASED"),
                age_d = ifelse(age < 65, 0, 1),
                stage = 0,
                nodes = 0,
                tumor = 0,
                metastasis = ifelse(ajcc_metastasis_pathologic_pm == "M0", 0, 1),
                smoke = "_yes") 

lusc$stage[grep("III|IV", lusc$ajcc_pathologic_tumor_stage)] <- 1

lusc$nodes[grep("2|3", lusc$ajcc_nodes_pathologic_pn)] <- 1

lusc$tumor[grep("3|4", lusc$ajcc_tumor_pathologic_pt)] <- 1

lusc$smoke[grep("3|4|5", lusc$tobacco_smoking_history_indicator)] <- "_reformed"
lusc$smoke[grep("1", lusc$tobacco_smoking_history_indicator)] <- "_no"

# An International Prognostic Index variable is created which has levels low, medium and high. 
lusc <- lusc %>%
  dplyr::mutate(ipi_n = stage + nodes + 
                  tumor + metastasis,
                ipi = ifelse( ipi_n < 1 , "low", 
                              ifelse( ipi_n < 2 , "medium",
                                      "hihg") ) )

lusc <- lusc %>%
  dplyr::mutate_at(c("stage", "nodes", "tumor", "smoke"), as.factor) %>%
  dplyr::select("patient_id","os_months", "os_status", "age", "stage", "nodes", "tumor", "metastasis", "smoke", "ipi")

lusc <- lusc[complete.cases(lusc), ] 

preProcValues <- caret::preProcess(lusc[c("age")], method = c("center", "scale") )
trainTransformed <- as.data.frame( predict(preProcValues, lusc[c("age")]) )
colnames(trainTransformed) = c("age_std")
luscs <- cbind(lusc, trainTransformed)
rm(list = c("trainTransformed", "preProcValues")) # but we keep preProcValues

#Train null model

system.time( survbayes.mod1 <- post_surv(x = luscs %>%
                                           dplyr::mutate(time = os_months,
                                                         status = os_status)) )

#To create step stan function
deparse(survbayes.mod1$formula)
survbayes.mod1$stan_function
#Stepwise selection
system.time( survbayes.step.fit1 <- stan_surv_step(x = luscs %>%
                                                     dplyr::mutate(time = os_months,
                                                                   status = os_status), 
                                                   fit = survbayes.mod1,
                                                   scope = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis", "s(age_std)"), verbose = TRUE)
)

system.time( survbayes.step.fit2 <- stan_surv_step(x = luscs %>%
                                                     dplyr::mutate(time = os_months,
                                                                   status = os_status), 
                                                   fit = survbayes.mod1,
                                                   scope = c("smoke", "age_std", "ipi", "s(age_std)"), verbose = TRUE)
)

pem.fits <- list(survbayes.step.fit1, survbayes.step.fit2)
pem.fits.dev <- sapply(fits, function(m){
  -2*m[[2]]$estimates["elpd_loo","Estimate"]
});max(unlist(fits.dev))

best.pem.fit <- pem.fits[[ match( max(unlist(pem.fits.dev)), pem.fits.dev ) ]][[1]]
class(best.pem.fit)

## Train Cox model
cox.mod1 <- cox_step(x = luscs %>% 
                       mutate(time = os_months,
                              status = os_status), surv_form = c("smoke", "age_std", "stage", "nodes", "tumor", "metastasis") )

cox.mod2 <- cox_step(x = luscs %>% 
                       mutate(time = os_months,
                              status = os_status), surv_form = c("smoke", "age_std", "ipi", "s(age_std)" ) )

cox.fits <- list(cox.mod1, cox.mod2)
cox.aic <- sapply(cox.fits , function(c){
  -2*c$loglik[[2]] + length( c$coefficients )
})

best.cox.fit <- cox.fits[[ match( max(unlist(cox.aic)), cox.aic ) ]]
class(best.cox.fit)

# Train test validation
# set.seed(99)
# slusc <- resample_stratified(lusc, strata = c( "age_d", "smoke", "ipi"), sizeto = 100)

set.seed(10)
vfold <- rsample::mc_cv(lusc, times = 1, strata = "os_status")$splits$`1`
train_lusc <- rsample::analysis(vfold)
test_lusc <- rsample::assessment(vfold)

preProcValues <- caret::preProcess(train_lusc[c("age")], method = c("center", "scale") )
trainTransformed <- as.data.frame( predict(preProcValues, train_lusc[c("age")]) )
colnames(trainTransformed) = c("age_std")
train_luscs <- cbind(train_lusc, trainTransformed)
rm(trainTransformed) # but we keep preProcValues
testTransformed <- as.data.frame( predict(preProcValues, test_lusc[c("age")]) )
colnames(testTransformed) = c("age_std")
test_luscs <- cbind(test_lusc, testTransformed)
rm(list = c("preProcValues", "testTransformed")) 

pem.train <- post_surv(x = train_luscs %>%
                         dplyr::mutate(time = os_months,
                                       status = os_status),
                       form = best.pem.fit$formula) 

#Train cox model
cox.train <- coxph( best.cox.fit$formula, data = train_luscs %>%
                      dplyr::mutate(time = os_months,
                                    status = os_status) )



## Schemper and Henderson Index
survAUC::schemper(train.fit = cph(Surv(time, status) ~ tumor + metastasis +  smoke + age_std, data = train_luscs %>%
                                    dplyr::mutate(time = os_months,
                                                  status = os_status)), traindata = train_luscs %>%
                    dplyr::mutate(time = os_months,
                                  status = os_status),
                  newdata = test_luscs %>%
                    dplyr::mutate(time = os_months,
                                  status = os_status))


## Royston D index
risk.prediction = predict(cox.train, newdata = test_luscs %>%
                            dplyr::mutate(time = os_months,
                                          status = os_status),
                          predict = "lp")

require(bootstrap)
require(KernSmooth)
survcomp::D.index(x = risk.prediction, surv.time = test_luscs$os_months, surv.event = test_luscs$os_status, strat = )

## Obtain Brier score

cox.brier <- get_survbrier(test = test_luscs %>%
                             dplyr::mutate(time = os_months,
                                           status = os_status),
                           mod = cox.train
)
cox.apperror <- cox.brier$AppErr$matrix

pem.brier <- get_survbrier(test = test_luscs %>%
                             dplyr::mutate(time = os_months,
                                           status = os_status),
                           mod = pem.train
)


## Obtain c-index

cox.cindex <- get_cindex(test = test_luscs %>%
                           dplyr::mutate(time = os_months,
                                         status = os_status),
                         mod = cox.train
)

brier.cindex <- get_cindex(test = test_luscs %>%
                             dplyr::mutate(time = os_months,
                                           status = os_status),
                           mod = pem.train
)


## Obtain ROC


cox.roc <- get_survroc(test = test_luscs %>%
                         dplyr::mutate(time = os_months,
                                       status = os_status),
                       mod = cox.train
)

pem.roc <- get_survroc(test = test_luscs %>%
                         dplyr::mutate(time = os_months,
                                       status = os_status),
                       mod = pem.train
)

end_time <- Sys.time()
end_time - start_time

