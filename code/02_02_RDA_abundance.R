## ------------------------ Script for analyzing Taxa abundances with RDA ------------------- ##

###### Section 1 ######
library(dplyr)
library(ggplot2)
library(vegan)
library(reshape2)

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

###### Section 2 Load abundance data ######

##/////////////////////////////////
## ------- June data ------------
##/////////////////////////////////
## abundances on family level
fam_jn <- read.csv(file.path(datadir, "invertebrate_data/familyjune_wide.csv"), stringsAsFactors = FALSE)
rownames(fam_jn) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "K", 
                      "L", "M", "N", "O", "P", "Q", "R", "S", "T")

names(fam_jn) <- gsub("^(\\w*)(.*?)$", "\\1", names(fam_jn))
## check gradient length
decorana(fam_jn)
## gradient lenght is fine

###### Section 3 Load environmental data ######
env_jn <- read.csv(file.path(datadir, "environmental_data/env_predictors_june.csv"), stringsAsFactors = FALSE)
## drop site column
env_jn <- env_jn[,-c(1)]

## run pca on environmental data
pca_env <- rda(env_jn, scale = TRUE)
summary(pca_env)
plot(pca_env, scaling = 1)

## check relationships between most abundant taxa and environmental variables
cor_dat <- cbind(env_jn, fam_jn[,c(2,6,13,20,27,31)])
pairs(cor_dat[,c(1:11)], lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(cor_dat[,c(1:10,12)], lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(cor_dat[,c(1:10,13)], lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(cor_dat[,c(1:10,14)], lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(cor_dat[,c(1:10,15)], lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(cor_dat[,c(1:10,16)], lower.panel = panel.smooth, upper.panel = panel.cor)
pairs(cor_dat[,c(1:10,17)], lower.panel = panel.smooth, upper.panel = panel.cor)


###### Section 4 RDA Model building ######
jn_rda0 <- rda(fam_jn ~ 1, data = env_jn, scale = TRUE)
jn_rda1 <- rda(fam_jn ~ ., data = env_jn, scale = TRUE)
## Forward selection using P-value as selection method
famtax_step <- ordiR2step(jn_rda0, scope = formula(jn_rda1), direction = "both", pstep = 1000, Pin = 0.1, R2scope = FALSE)

## most parsimonious model
jfamtax_rda_pars <- rda(fam_jn ~ arableLand + pastures, data = env_jn, scale = TRUE)
summary(jfamtax_rda_pars)
## constrained model explains 18.7% of variation
## first two axes explain 18.7% of variation
scores(jfamtax_rda_pars, c(1,2))

## check model
set.seed(187)
sig_june <- (anova.cca(jfamtax_rda_pars, step = 1000))
## model is significant

# check axes of model
set.seed(187)
anova.cca(jfamtax_rda_pars, step = 1000, by = "axis")
## only first axis is barely significant

set.seed(187)
anova.cca(jfamtax_rda_pars, step = 1000, by = "term")
## first two terms are significant, third one not


vif.cca(jfamtax_rda_pars)
## Plot model with scaling 1: distance triplot
plot(jfamtax_rda_pars, scaling = 1, main = "Family Taxa June, RDA scaling 1 - wa scores")
# famtax_pars_sc <- scores(famtax_rda_pars, choices = 1:2, scaling = 1, display = "sp")
# arrows(0, 0, famtax_pars_sc[, 1], famtax_pars_sc[, 2], length = 0, lty = 1, col = "red")
s2 <- summary(jfamtax_rda_pars)


###### Section 5 Write RDA results in to a table ######
rda_tab <- data.frame(var = as.character(c("total inertia", "constrained variance", "unconstrained variance",
                                           "RDA 1 % explained", "RDA 2 % explained", "p", "observations")),
                      RDA_S2 = as.numeric(c(s2$tot.chi, s2$constr.chi/s2$tot.chi*100, s2$unconst.chi/s2$tot.chi*100, s2$cont[[1]][2,1]*100,
                                            s2$cont[[1]][2,2]*100, sig_june$`Pr(>F)`[1], nrow(s2$sites))),
                      stringsAsFactors = FALSE)
## round the values 
rda_tab <- rda_tab %>% 
  mutate(across(where(is.numeric), round, 2))
#rda_tab$RDA_S1 <- as.character(rda_tab$RDA_S1)
rda_tab$RDA_S2 <- as.character(rda_tab$RDA_S2)

## write table to file
write.csv(rda_tab, "./model_outputs/rda_tab.csv", row.names = FALSE)

###### Section 6 Plotting the results with ggplot ######
library(ggrepel)

summary(jfamtax_rda_pars)

ii_jn <- summary(jfamtax_rda_pars)
sp2 <- data.frame(scores(jfamtax_rda_pars, choices = c(1,2), display = c("species"), scaling = 1), stringsAsFactors = FALSE)
sp2$Label <- row.names(sp2)
st2 <- data.frame(scores(jfamtax_rda_pars, choices = c(1,2), display = c("sites"), scaling = 1), stringsAsFactors = FALSE)
st2$Label <- row.names(st2)
yz2 <- data.frame(scores(jfamtax_rda_pars, choices = c(1,2), display = c("bp"), scaling = 1), stringsAsFactors = FALSE)
yz2$Label <- row.names(yz2)
const2 <- 4.973795
yz2[,1:2] <- yz2[,1:2] * const2

# grp_jn <- as.data.frame(c(rep("a", 4), rep("b", 4), rep("c", 4), rep("d", 4)))
# colnames(grp) = "group"
pars_rda_jn <- ggplot() +
  geom_hline(yintercept = 0, linetype = 1, size = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.5, alpha = 0.5) +
  geom_point(data = st2, aes(RDA1, RDA2), shape = 2, color = "grey40" ) + 
  geom_text_repel(data = st2, aes(RDA1, RDA2, label = row.names(st2)), color = "grey40", size = 2.5) +
  # geom_segment(data = sp2, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
  #              arrow = arrow(angle = 22.5, length = unit(0.35, "cm"),
  #              type = "closed"), linetype = 1, size = 0.6, colour = "red") +
  geom_point(data = sp2, aes(RDA1, RDA2), shape = 0, colour = "black") + 
  geom_text_repel(data = sp2, aes(RDA1, RDA2, label = row.names(sp2)), color = "black", size = 3) +
  geom_segment(data = yz2, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.15, "cm"),
                             type = "closed"), linetype = 1, size = 0.6, colour = "orangered3") +
  geom_text(data = yz2[1,], aes(RDA1, RDA2, label = "arable land"), color = "orangered3", size = 5,
            nudge_y = -0.1) +
  geom_text(data = yz2[2,], aes(RDA1, RDA2, label = row.names(yz2)[2]), color = "orangered3", size = 5,
            nudge_y = 0.1) +
  labs(x=paste("RDA 1 (explained variance: ", format(100 * ii_jn$cont[[1]][2,1], digits = 3), " %)", sep= ""),
       y=paste("RDA 2 (explained variance: ", format(100 * ii_jn$cont[[1]][2,2], digits = 2), " %)", sep = "")) +
  # guides(shape = guide_legend(title = NULL, colour = "black"),
  #  fill = guide_legend(title = NULL)) +
  theme_bw() #+ theme(panel.grid.minor = element_blank())

pars_rda_jn

ggsave(file= "./plots/Biplot_RDA_June.png", plot=pars_rda_jn, width = 7, height = 4.5)
ggsave(file= "./plots/Biplot_RDA_June.svg", plot=pars_rda_jn, width = 7, height = 4.5)
ggsave(file= "./plots/Biplot_RDA_June.tiff", plot=pars_rda_jn, width = 7, height = 4.5)
