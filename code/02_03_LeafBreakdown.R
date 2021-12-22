## ------------------- Script for analyzing Leaf Litter Breakdown --------------------- ##


###### Section 1 Loading packages and file paths ######
library(data.table)
library(dplyr)
library(readxl)

## load connections
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

###### Section 2 Load data for leaf breakdown, calculate decomposition rate k ######
## load data of coarse and fine mesh bags
mloss_fine <- read.csv("./data/leaf_breakdown/lld_fine.csv", stringsAsFactors = FALSE)
mloss_coarse <- read.csv("./data/leaf_breakdown/lld_coarse.csv", stringsAsFactors = FALSE)
head(mloss_fine)
head(mloss_coarse)
## calculate mass loss in g for fine leaf bags
mloss_fine$mlf_g <- mloss_fine$mean_mass_start - mloss_fine$mean_mass_end
## merge mlf_g to mloss_coarse
mloss_coarse <- merge(mloss_coarse, mloss_fine[,c(1,2,8)], by = c("site", "leaf_type"))


## load environmental data to get temperatures
data <- read_xlsx(file.path(datadir, "environmental_data/raw_env_dat.xlsx"), sheet = 8)
## keep only temperature
data <- data[data$Variable == "water temperature" & data$Date >= "2016-05-22", c(3,4,6)]
## change Value column to numeric
data$Value <- as.numeric(data$Value)
## calculate mean temperatures
mean_temp <- data %>%
  group_by(Site) %>%
  summarize(mean_temp = mean(Value, na.rm = TRUE))
names(mean_temp)[1] <- "site"

## merge mean_temp to mloss_coarse
mloss_coarse <- merge(mloss_coarse, mean_temp, by = "site", all.x = TRUE)

## calculate kshred
mloss_coarse$kshred <- (-log((mloss_coarse$mean_mass_end + mloss_coarse$mlf_g)/mloss_coarse$mean_mass_start))/mloss_coarse$mean_temp
## cast into wide format
kshred_wide <- dcast(mloss_coarse[mloss_coarse$leaf_type %in% c("E", "B"),c(1,2,11)],
                     site ~ leaf_type, value.var = "kshred")
## replace negative values with 0
kshred_wide[c(9,15),2] <- c(0,0)
kshred_wide[c(17),3] <- 0
names(kshred_wide) <- c("site", "kb", "k")
# write.csv(kshred_wide[,c(1,3)], "./data/leaf_breakdown/lld_rate.csv", row.names = FALSE)

## extra set: Correlate SPEAR with k
## load macroinv_metrics
metrics <- read.csv("./data/invertebrate_data/macroinv_metrics.csv", stringsAsFactors = FALSE)
## subset to SPEAR
metrics <- metrics[metrics$month == "June", c(1,12)]
## merge with kshred
metrics <- merge(metrics, kshred_wide, by = "site")
plot(metrics$spear_pest ~ metrics$k)
cor(metrics$spear_pest, metrics$k)
cor.test(metrics$spear_pest, metrics$k)

###### Section 3 Load environmental data and shredder data ######
## load environmental data
env_jn <- read.csv(file.path(datadir, "environmental_data/env_predictors_june.csv"), stringsAsFactors = FALSE)

## load Gammarid-biomass
shr_com <- read.csv(file.path(datadir, "invertebrate_data/shredder_community.csv"), stringsAsFactors = FALSE)
shr_com <- shr_com[,-c(6:8)]
## round gam_bodyweight
shr_com$gam_bodyweight <- round(shr_com$gam_bodyweight, 3)
shr_com$log_gamshr <- log10(shr_com$gam_shr + 1)
shr_com$log_allshr <- log10(shr_com$all_shr)
shr_com$log_nongam <- log10(shr_com$nongam_shr + 1)

## merge all dataframes
shred_dat <- merge(shr_com[,-c(6:8)], env_jn, by = "site", all = TRUE)
massl_dat <- merge(kshred_wide[,c(1,3)], env_jn, by = "site", all = TRUE)
massl_dat <- merge(massl_dat, shr_com[,-c(3:5)], by = "site", all = TRUE)
## drop site M
massl_dat <- massl_dat[-c(12),]


###### Section 4 Check datasets ######
## use pairs plot
## shredder against environmental drivers
pairs(shred_dat[,-1], lower.panel = panel.smooth, upper.panel = panel.cor)
shred_dat$log_gam <- log10(shred_dat$gam_shr + 1)
shred_dat$log_shr <- log10(shred_dat$all_shr + 1)
shred_dat$log_ng <- log10(shred_dat$nongam_shr + 1)
pairs(shred_dat[,c(16:18,6:15)], lower.panel = panel.smooth, upper.panel = panel.cor)

## mass loss against environmental drivers and shredders
pairs(massl_dat[,-1], lower.panel = panel.smooth, upper.panel = panel.cor)



###### Section 5 Model shredder abundance and body weight of Gammarus balcanicus ######
# library(MASS)
library(glmnet)
library(caret)


## Aim of elastic net regression is to identify the most relevant variables
## Lasso regression shrinks regression coefficients of unimportant variables to zero

##///////////////////////////////////////////////////
## ------ 5.1 body weight as dependent variable -----
##///////////////////////////////////////////////////
## subset data
bw_data <- shred_dat[,c(2,6:15)]
## remove 0-entry of bodyweight, is positive continuous variable
bw_data <- bw_data[-11,]

## /////////// continue with elastic net regression and caret package
# standardise exp. variables
x <- bw_data %>% select(-gam_bodyweight) %>% as.matrix()
y <- bw_data %>% select(gam_bodyweight) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
bowe_elnet_mod <- train(gam_bodyweight ~ .,
                         data = cbind(y, x),
                         method = "glmnet",
                         family = "gaussian",
                         metric = "RMSE",
                         preProcess = c("center", "scale"),
                         tuneLength = 50,
                         trControl = train_control)

## function to extract row with best tuning parameters
res_bowe_elnet <- data.frame(get_best_result(bowe_elnet_mod))
res_bowe_elnet$mod <- "Avergae body weight"
coef_bowe_elnet <- data.frame(as.matrix(coef(bowe_elnet_mod$finalModel, bowe_elnet_mod$bestTune$lambda)))
## plot final model
plot(bowe_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(bowe_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## standardize variables
x1 <- bw_data %>% select(-gam_bodyweight) %>% scale() %>% as.matrix()
y1 <- bw_data %>% select(gam_bodyweight) %>% as.matrix()
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
bowe_mod <- glmnet(x1, y1, alpha = res_bowe_elnet$alpha, lambda = res_bowe_elnet$lambda, family = gaussian())
coef.glmnet(bowe_mod)
## get regression coefficients for model with best lambda
coef_bowe_elnet <- data.frame(as.matrix(coef.glmnet(bowe_mod)))
names(coef_bowe_elnet) <- "bowe"
bowe_r2 <- bowe_mod$dev.ratio
bowe_lam <- res_bowe_elnet$lambda
bowe_alpha <- res_bowe_elnet$alpha



##///////////////////////////////////////////////////
## ------ 5.2 Abundance of Gammarus balcanicus  -----
##///////////////////////////////////////////////////

## overall shredder abundance is dominated by Gammarids
plot(shred_dat$all_shr ~ shred_dat$gam_shr)
## results from all_shred also apply for Gammarid shredders

## subset data
ab_data <- shred_dat[,c(16,6:15)]

## /////////// continue with elastic net regression and caret package
# standardize exp. variables
x <- ab_data %>% select(-log_gam) %>% as.matrix()
y <- ab_data %>% select(log_gam) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
gash_elnet_mod <- train(log_gam ~ .,
                        data = cbind(y, x),
                        method = "glmnet",
                        family = gaussian(),
                        metric = "RMSE",
                        preProcess = c("center", "scale"),
                        tuneLength = 50,
                        trControl = train_control)

## function to extract row with best tuning parameters
res_gash_elnet <- data.frame(get_best_result(gash_elnet_mod))
res_gash_elnet$mod <- "Gammarid abundance"
## plot final model
plot(gash_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(gash_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## standardize variables
x1 <- ab_data %>% select(-log_gam) %>% scale() %>% as.matrix()
y1 <- ab_data %>% select(log_gam) %>% as.matrix()
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
gash_mod <- glmnet(x1, y1, alpha = res_gash_elnet$alpha, lambda = res_gash_elnet$lambda, family = gaussian())
coef.glmnet(gash_mod)
## get regression coefficients for model with best lambda
coef_gash_elnet <- data.frame(as.matrix(coef.glmnet(gash_mod)))
names(coef_gash_elnet) <- "log_gash"
gash_r2 <- gash_mod$dev.ratio
gash_lam <- res_gash_elnet$lambda
gash_alpha <- res_gash_elnet$alpha


##///////////////////////////////////////////////////
## --- 5.3 Abundance of Non-Gammarid shredders  ----
##///////////////////////////////////////////////////

## subset data
ng_data <- shred_dat[,c(18,6:15)]

## /////////// continue with elastic net regression and caret package
# standardize exp. variables
x <- ng_data %>% select(-log_ng) %>% as.matrix()
y <- ng_data %>% select(log_ng) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
shr_elnet_mod <- train(log_ng ~ .,
                        data = cbind(y, x),
                        method = "glmnet",
                        family = gaussian(),
                        metric = "RMSE",
                        preProcess = c("center", "scale"),
                        tuneLength = 50,
                        trControl = train_control)

## function to extract row with best tuning parameters
res_shr_elnet <- data.frame(get_best_result(shr_elnet_mod))
res_shr_elnet$mod <- "log non-Gammarid Shredder abundance"
## plot final model
plot(shr_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(shr_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## standardize variables
x1 <- ng_data %>% select(-log_ng) %>% scale() %>% as.matrix()
y1 <- ng_data %>% select(log_ng) %>% as.matrix()
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
shr_mod <- glmnet(x1, y1, alpha = res_shr_elnet$alpha, lambda = res_shr_elnet$lambda, family = gaussian())
coef.glmnet(shr_mod)
## get regression coefficients for model with best lambda
coef_shr_elnet <- data.frame(as.matrix(coef.glmnet(shr_mod)))
names(coef_shr_elnet) <- "shr"
shr_r2 <- shr_mod$dev.ratio
shr_lam <- shr_mod$lambda
shr_alpha <- shr_mod$alpha


###### Section 6 Model leaf litter breakdown ######
##///////////////////////////////////////////////////
## ------ 6.1 k as response variable  --------------
##///////////////////////////////////////////////////
k_data <- massl_dat[,c(2:12,14,16)]

## /////////// continue with elastic net regression and caret package
# standardise exp. variables
x <- k_data %>% select(-k) %>% as.matrix()
y <- k_data %>% select(k) %>% as.matrix()

## Set parameters for training control
train_control <- trainControl(method = "LOOCV", ## small number of observations, use LOOCV
                              search = "grid",
                              verboseIter = TRUE)

## Train the model
k_elnet_mod <- train(k ~ .,
                        data = cbind(y, x),
                        method = "glmnet",
                        family = "gaussian",
                        metric = "RMSE",
                        preProcess = c("center", "scale"),
                        tuneLength = 50,
                        trControl = train_control)

## function to extract row with best tuning parameters
res_k_elnet <- data.frame(get_best_result(k_elnet_mod))
res_k_elnet$mod <- "k"
plot(k_elnet_mod$finalModel, xvar = "lambda", label = TRUE)
plot(k_elnet_mod$finalModel, xvar = "dev", label = TRUE)

## //////// refit the model, get the deviance ratio (r2)
## standardize variables
x1 <- k_data %>% select(-k) %>% scale() %>% as.matrix()
y1 <- k_data %>% select(k) %>% as.matrix()
## according to https://stackoverflow.com/questions/50610895/how-to-calculate-r-squared-value-for-lasso-regression-using-glmnet-in-r
k_mod <- glmnet(x1, y1, alpha = res_k_elnet$alpha, lambda = res_k_elnet$lambda, family = gaussian())
coef.glmnet(k_mod)
## get regression coefficients for model with best lambda
coef_k_elnet <- data.frame(as.matrix(coef.glmnet(k_mod)))
names(coef_k_elnet) <- "k"
k_r2 <- k_mod$dev.ratio
k_lam <- res_k_elnet$lambda
k_alpha <- res_k_elnet$alpha


###### Section 7 Combine results of elastic net models ######

## ===>> Combine results of elastic net models
model_parameters <- data.frame(Variable = as.character(c("alpha", "lambda", "r2")),
                               bowe = as.numeric(c(bowe_alpha, bowe_lam, bowe_r2)),
                               gash = as.numeric(c(gash_alpha, gash_lam, gash_r2)),
                               nong = as.numeric(c(0, shr_lam, shr_r2)),
                               k = as.numeric(c(k_alpha, k_lam, k_r2)),
                               stringsAsFactors = FALSE)
names(model_parameters) <- c("Variable", "bowe", "log_gash", "log_nongam", "k")
model_parameters <- model_parameters %>% mutate(across(where(is.numeric), round, 3))
#write.csv(model_parameters, "./model_outputs/model_para.csv", row.names = FALSE)

## ===>> Combine coefficients of elastic net models
df1 <- data.frame(bowe = as.numeric(c(NA,NA)),
                  row.names = c("log_gamshr", "log_nongam"), 
                  stringsAsFactors = FALSE)
df2 <- data.frame(log_gash = as.numeric(c(NA, NA)),
                  row.names = c("log_gamshr", "log_nongam"), 
                  stringsAsFactors = FALSE)
df3 <- data.frame(shr = as.numeric(c(NA, NA)),
                  row.names = c("log_gamshr", "log_nongam"),
                  stringsAsFactors = FALSE)
coef_bowe_elnet <- rbind(coef_bowe_elnet, df1)
coef_gash_elnet <- rbind(coef_gash_elnet, df2)
coef_shr_elnet <- rbind(coef_shr_elnet, df3)
model_coefficients <- cbind(coef_bowe_elnet, coef_gash_elnet, coef_shr_elnet, coef_k_elnet)
model_coefficients$Variable <- row.names(model_coefficients)
model_coefficients <- model_coefficients[,c(5,1:4)]
names(model_coefficients) <- c("Variable", "bowe", "log_gash", "log_nongam", "k")
model_coefficients <- model_coefficients %>% mutate(across(where(is.numeric), round, 2))
#write.csv(model_coefficients, "./model_outputs/model_coef.csv", row.names = FALSE)

## combine both tables
model_table <- rbind(model_coefficients, model_parameters)
model_table[] <- lapply(model_table, as.character)
family_df <- data.frame(Variable = "family",
                        bowe = "Gaussian",
                        log_gash = "Gaussian",
                        log_nongam = "Gaussian",
                        k = "Gaussian")
model_table <- rbind(model_table, family_df)
## change order of columns
model_table <- model_table[c(17,16,14,15,1,9,8,11,10,7,4,6,2,5,3,12,13),]
write.csv(model_table, "./model_outputs/model_shred_tab.csv", row.names = FALSE)

