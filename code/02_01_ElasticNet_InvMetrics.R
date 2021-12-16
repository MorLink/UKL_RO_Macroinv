## -------------------------------- Script for calculating SPEAR values -------------------------- ##

###### Section 1 Loading packages and file paths ######
## --------------- Load packages -------------------##
library(data.table)
library(dplyr)
library(ggplot2)
library(glmnet)
library(plotmo)
library(stabs)
library(effects)
library(caret)

## -----------------load paths ----------------------##
wd <- setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'


## ---------------- load functions ---------------##
## check correlation plots
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex <- 1 / strwidth(txt)
  
  test <- cor.test(x, y)
  # borrowed from printCoefmat
  Signif <- symnum(
    test$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.05, 0.1, 1),
    symbols = c("*", ".", " ")
  )
  
  text(0.5, 0.5, txt, cex = cex * r)
  text(5, 5, Signif, cex = cex, col = 2)
}

## function to extract row with best tuning parameters
get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

###### Section 2 Load macroinvertebrate metrics and environmental data ######
##//////////////////////////////////////////
## ----- Load macroinvertebrate metrics ----
##//////////////////////////////////////////
inv_metrics <- read.csv(file.path(datadir, "invertebrate_data/macroinv_metrics.csv"), stringsAsFactors = FALSE)

##//////////////////////////////////////////
## ----- Load environmental data ----
##//////////////////////////////////////////
env_jn <- read.csv(file.path(datadir, "environmental_data/env_predictors_june.csv"), stringsAsFactors = FALSE)

## merge invertebrate metrics and environmental data
dat_jn <- merge(env_jn, inv_metrics[inv_metrics$month == "June",-2], by = "site")


###### Section 3 Modeling the invertebrate metrics ######

##/////////////////////////////////////////////
## ---- 3.1 SWD as dependent variable ---------
##/////////////////////////////////////////////
swd_cols <- c(2:11,14)
swd_dat <- dat_jn[swd_cols]

## check pairs plot
pairs(swd_dat, lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()

## ///// use stepwise forward model building with BIC
reg0 <- lm(swd ~ 1,
           data = swd_dat)
reg1 <- lm(swd ~ .,
           data = swd_dat)
step(reg0, scope=formula(reg0, reg1),
     direction="forward", k = log(nrow(swd_dat)))
## -->> no independent variable selected

## standardize variables
x1 <- swd_dat %>% select(-swd) %>% scale() %>% as.matrix()
y1 <- swd_dat %>% select(swd) %>% as.matrix()

## ///// use cv.glmnet for feature selection
set.seed(4159)
cv_tnr <- cv.glmnet(x1, y1, alpha = 1, nfolds = 5, family = "gaussian")  ## 19 obs, 10 folds not possible, 3 folds too small
plot(cv_tnr)
coef(cv_tnr, s = "lambda.1se") ## all coefficients shrunk to zero
coef(cv_tnr, s = "lambda.min") ## fine_sed selected

## ///// Apply stability selection with stabs package
set.seed(3254)
(stabs.elnet <- stabsel(x = x, y = y, fitfun = glmnet.lasso,
                        cutoff = 0.60, PFER = 1))
## no variables selected

## /////////// continue with elastic net regression and caret package
# standardise exp. variables
x <- swd_dat %>% select(-swd) %>% as.matrix()
y <- swd_dat %>% select(swd) %>% scale(center = TRUE, scale = FALSE) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
swd_elnet_mod <- train(swd ~ .,
                           data = cbind(y, x),
                           method = "glmnet",
                           family = "gaussian",
                           preProcess = c("center", "scale"),
                           tuneLength = 50,
                           trControl = train_control)

## function to extract row with best tuning parameters
res_swd_elnet <- data.frame(get_best_result(swd_elnet_mod))
res_swd_elnet$mod <- "swd"
## plot final model
plot(swd_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(swd_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
swd_mod <- glmnet(x1, y1, alpha = res_swd_elnet$alpha, lambda = res_swd_elnet$lambda, family = gaussian())
coef.glmnet(swd_mod)
## get regression coefficients for model with best lambda
coef_swd_elnet <- data.frame(as.matrix(coef.glmnet(swd_mod)))
names(coef_swd_elnet) <- "swd"
swd_r2 <- swd_mod$dev.ratio
swd_lam <- res_swd_elnet$lambda
swd_alpha <- res_swd_elnet$alpha


##///////////////////////////////////////////////////
## - 3.2 Percent EPT taxa as dependent variable ----
##///////////////////////////////////////////////////
epttp_cols <- c(2:11,17)
epttp_dat <- dat_jn[epttp_cols]

## check pairs plot
pairs(epttp_dat, lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()

## ///// use stepwise forward model building with BIC
reg0 <- lm(epttax_percent ~ 1,
           data = epttp_dat)
reg1 <- lm(epttax_percent ~ .,
           data = epttp_dat)
step(reg0, scope=formula(reg0, reg1),
     direction="forward", k = log(nrow(epttp_dat)))
## -->> no independent variable selected

## standardize variables
x1 <- epttp_dat %>% select(-epttax_percent) %>% scale() %>% as.matrix()
y1 <- epttp_dat %>% select(epttax_percent) %>% as.matrix()

## ///// use cv.glmnet for feature selection
set.seed(4159)
cv_tnr <- cv.glmnet(x1, y1, alpha = 1, nfolds = 5, family = "gaussian")  ## 19 obs, 10 folds not possible, 3 folds too small
plot(cv_tnr)
coef(cv_tnr, s = "lambda.1se") ## all coefficients shrunk to zero
coef(cv_tnr, s = "lambda.min") ## water_temp, fine_sed, agriculture, refugium

## ///// Apply stability selection with stabs package
set.seed(3254)
(stabs.elnet <- stabsel(x = x1, y = y1, fitfun = glmnet.lasso,
                        cutoff = 0.60, PFER = 1))
## no variables selected 

## /////////// continue with elastic net regression and caret package
# standardise exp. variables
x <- epttp_dat %>% select(-epttax_percent) %>% as.matrix()
y <- epttp_dat %>% select(epttax_percent) %>% scale(center = TRUE, scale = FALSE) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
epttp_elnet_mod <- train(epttax_percent ~ .,
                        data = cbind(y, x),
                        method = "glmnet",
                        family = "gaussian",
                        preProcess = c("center", "scale"),
                        tuneLength = 50,
                        trControl = train_control)

## function to extract row with best tuning parameters
res_epttp_elnet <- data.frame(get_best_result(epttp_elnet_mod))
res_epttp_elnet$mod <- "epttax_percent"


## plot final model
plot(epttp_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(epttp_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
epttp_mod <- glmnet(x1, y1, alpha = res_epttp_elnet$alpha, lambda = res_epttp_elnet$lambda, family = gaussian())
coef.glmnet(epttp_mod)
## get regression coefficients for model with best lambda
coef_epttp_elnet <- data.frame(as.matrix(coef.glmnet(epttp_mod)))
names(coef_epttp_elnet) <- "epttax_percent"
epttp_r2 <- epttp_mod$dev.ratio
epttp_lam <- res_epttp_elnet$lambda
epttp_alpha <- res_epttp_elnet$alpha


##///////////////////////////////////////////////////
## - 3.3 Percent EPT abundance as dependent variable ----
##///////////////////////////////////////////////////
eptp_cols <- c(2:11,18)
eptp_dat <- dat_jn[eptp_cols]

## check pairs plot
pairs(eptp_dat, lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()

## ///// use stepwise forward model building with BIC
reg0 <- lm(ept_percent ~ 1,
           data = eptp_dat)
reg1 <- lm(ept_percent ~ .,
           data = eptp_dat)
step(reg0, scope=formula(reg0, reg1),
     direction="forward", k = log(nrow(eptp_dat)))
## -->> no independent variable selected

## standardize variables
x1 <- eptp_dat %>% select(-ept_percent) %>% scale() %>% as.matrix()
y1 <- eptp_dat %>% select(ept_percent) %>% as.matrix()

## ///// use cv.glmnet for feature selection
set.seed(4159)
cv_tnr <- cv.glmnet(x1, y1, alpha = 1, nfolds = 5, family = "gaussian")  ## 19 obs, 10 folds not possible, 3 folds too small
plot(cv_tnr)
coef(cv_tnr, s = "lambda.1se") ## all coefficients shrunk to zero
coef(cv_tnr, s = "lambda.min") ## water_temp, fine_sed, agriculture, refugium

## ///// Apply stability selection with stabs package
set.seed(3254)
(stabs.elnet <- stabsel(x = x, y = y, fitfun = glmnet.lasso,
                        cutoff = 0.60, PFER = 1))
## no variables selected

## /////////// continue with elastic net regression and caret package
# standardise exp. variables
x <- eptp_dat %>% select(-ept_percent) %>% as.matrix()
y <- eptp_dat %>% select(ept_percent) %>% scale(center = TRUE, scale = FALSE) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
eptp_elnet_mod <- train(ept_percent ~ .,
                           data = cbind(y, x),
                           method = "glmnet",
                           family = "gaussian",
                           preProcess = c("center", "scale"),
                           tuneLength = 50,
                           trControl = train_control)

## function to extract row with best tuning parameters
res_eptp_elnet <- data.frame(get_best_result(eptp_elnet_mod))
res_eptp_elnet$mod <- "ept_percent"


## plot final model
plot(eptp_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(eptp_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
eptp_mod <- glmnet(x1, y1, alpha = res_eptp_elnet$alpha, lambda = res_eptp_elnet$lambda, family = gaussian())
coef.glmnet(eptp_mod)
## get regression coefficients for model with best lambda
coef_eptp_elnet <- data.frame(as.matrix(coef.glmnet(eptp_mod)))
names(coef_eptp_elnet) <- "ept_percent"
eptp_r2 <- eptp_mod$dev.ratio
eptp_lam <- res_eptp_elnet$lambda
eptp_alpha <- res_eptp_elnet$alpha


##/////////////////////////////////////////////
## - 3.4 SPEARpest as dependent variable -----
##/////////////////////////////////////////////
spear_cols <- c(2:11,21)
spear_dat <- dat_jn[spear_cols]

## check pairs plot
pairs(spear_dat, lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()

## ///// use stepwise forward model building with BIC
reg0 <- lm(spear_pest ~ 1,
           data = spear_dat)
reg1 <- lm(spear_pest ~ .,
           data = spear_dat)
step(reg0, scope=formula(reg0, reg1),
     direction="forward", k = log(nrow(spear_dat)))
## -->> no independent variable selected

## standardize variables
x1 <- spear_dat %>% select(-spear_pest) %>% scale() %>% as.matrix()
y1 <- spear_dat %>% select(spear_pest) %>% as.matrix()

## ///// use cv.glmnet for feature selection
set.seed(4159)
cv_tnr <- cv.glmnet(x1, y1, alpha = 1, nfolds = 5, family = "gaussian")  ## 19 obs, 10 folds not possible, 3 folds too small
plot(cv_tnr)
coef(cv_tnr, s = "lambda.1se") ## all coefficients shrunk to zero
coef(cv_tnr, s = "lambda.min") ## refugium, flow, max_sumtu_iv


## ///// Apply stability selection with stabs package
set.seed(3254)
(stabs.elnet <- stabsel(x = x, y = y, fitfun = glmnet.lasso,
                        cutoff = 0.60, PFER = 1))
## no variables selected

## /////////// continue with elastic net regression and caret package
# standardise exp. variables
x <- spear_dat %>% select(-spear_pest) %>% as.matrix()
y <- spear_dat %>% select(spear_pest) %>% scale(center = TRUE, scale = FALSE) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
spear_elnet_mod <- train(spear_pest ~ .,
                           data = cbind(y, x),
                           method = "glmnet",
                           family = "gaussian",
                           metric = "RMSE",
                           preProcess = c("center", "scale"),
                           tuneLength = 50,
                           trControl = train_control)

## function to extract row with best tuning parameters
res_spear_elnet <- data.frame(get_best_result(spear_elnet_mod))
res_spear_elnet$mod <- "spear_pest"
## plot final model
plot(spear_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(spear_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
spear_mod <- glmnet(x1, y1, alpha = res_spear_elnet$alpha, lambda = res_spear_elnet$lambda, family = gaussian())
coef.glmnet(spear_mod)
## get regression coefficients for model with best lambda
coef_spear_elnet <- data.frame(as.matrix(coef.glmnet(spear_mod)))
names(eptp_spear_elnet) <- "spear_pest"
spear_r2 <- spear_mod$dev.ratio
spear_lam <- res_spear_elnet$lambda
spear_alpha <- res_spear_elnet$alpha


###### Section 4: Combine results ######

## ===>> Combine results of elastic net models
model_parameters <- data.frame(Variable = as.character(c("alpha", "lambda", "r2")),
                               SWD = as.numeric(c(swd_alpha, swd_lam, swd_r2)),
                               Percent_EPT_Taxa = as.numeric(c(epttp_alpha, epttp_lam, epttp_r2)),
                               Percent_EPT = as.numeric(c(eptp_alpha, eptp_lam, eptp_r2)),
                               spear = as.numeric(c(spear_alpha, spear_lam, spear_r2)),
                               stringsAsFactors = FALSE)
names(model_parameters) <- c("Variable", "SWD", "Percent_EPT_Taxa", "Percent_EPT", "SPEAR")
model_parameters <- model_parameters %>% mutate(across(where(is.numeric), round, 2))
#write.csv(model_parameters, "./model_outputs/model_para.csv", row.names = FALSE)

## ===>> Combine coefficients of elastic net models
model_coefficients <- cbind(coef_swd_elnet, coef_epttp_elnet, coef_eptp_elnet, coef_spear_elnet)
model_coefficients$Variable <- row.names(model_coefficients)
model_coefficients <- model_coefficients[,c(5,1:4)]
names(model_coefficients) <- c("Variable", "SWD", "Percent_EPT_Taxa", "Percent_EPT", "SPEAR")
model_coefficients <- model_coefficients %>% mutate(across(where(is.numeric), round, 2))
# model_coefficients[] <- lapply(model_coefficients, as.character)
#write.csv(model_coefficients, "./model_outputs/model_coef.csv", row.names = FALSE)

## combine both tables
model_table <- rbind(model_coefficients, model_parameters)
model_table[] <- lapply(model_table, as.character)
family_df <- data.frame(Variable = "family",
                        SWD = "Gaussian",
                        Percent_EPT_Taxa = "Gaussian",
                        Percent_EPT = "Gaussian",
                        SPEAR = "Gaussian",
                        stringsAsFactors = FALSE)
model_table <- rbind(model_table, family_df)
write.csv(model_table, "./model_outputs/model_tab.csv", row.names = FALSE)
