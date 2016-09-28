# load packages and functions
# source("diagnostic_fcns.r")
# library(lme4)
# library(car)
# library(lattice)
# library(arm)
library(ggplot2)

## Load dataset
df <- read.csv("datasets/dataset_5.csv")

head(df)
str(df)

df$Individual <- as.factor(df$Individual)

########################################################################
############ Issues to consider before fitting a model #################
########################################################################

#### Data distribution
# distribution.plot(df$Phenotype)


################
# COLLINEARITY #
################

GGally::ggpairs(
	df, columns = c("X1", "X2", "Phenotype"),
	diag  = list(continuous = "barDiag", colour = "grey"),
	lower = list(
		continuous = "smooth",
		mapping    = aes(color = Individual)
	)
)

# VIF: Variance inflation factor
# if VIF < threshold, no collinearity
# Threshold =  10 (Montgomery, D.C. & Peck, E.A. 1992. Wiley)
# Threshold =  3 (Zuur, A.F. et al. 2010. Methods in Ecology and Evolution)

mod  <- lm(Phenotype ~ X1 + X2, data = df)
car::vif(mod)

# Collinearity can be solved by dropping collinear covariates.
# Using VIF
# or (prehaps better) use common sens and biological knowledge


################
# MISSING DATA #
################

# df$missing <- ifelse(is.na(df$Phenotype), TRUE, FALSE)
# 
# Amelia::missmap(df)
# 
# # test if predictors of missing data are different from the not missing data
# t.test(X1~missing, df)
# boxplot(X1~missing, df)


# # Response variable
# par(mfrow = c(2, 2))
# # Response variable
# boxplot(df$Phenotype, ylab = "Phenotype", main = "Boxplot")
# dotchart(df$Phenotype,
# 				 ylab = "Order of observations",
# 				 xlab = "Phenotype", main = "Cleveland dotplot")
# boxplot(df$Phenotype~df$Individual, ylab = "X1", xlab = "Individual", main = "Boxplot")
# dotchart(df$Phenotype,
# 				 groups = factor(df$Individual),
# 				 ylab = "Individual", xlab = "Phenotype",
# 				 main = "Cleveland dotplot", pch = df$Individual)
# 
# 
# # Explanatory variable
# par(mfrow = c(2, 2))
# # Explanatory variable
# boxplot(df$X1, ylab = "X1", main = "Boxplot")
# dotchart(df$X1, ylab = "Order of observations", xlab = "X1", main = "Cleveland dotplot")
# boxplot(df$X1~df$Individual, ylab = "X1", xlab = "Individual", main = "Boxplot")
# dotchart(df$X1,
# 				 groups = factor(df$Individual),
# 				 ylab = "Individual", xlab = "X1",
# 				 main = "Cleveland dotplot", pch = df$Individual)


###################
# MODEL SELECTION #
###################

# Zuur, A. F. et al. (2009). Mixed Effects Models and Extensions in Ecology with R.


## Step 1: beyond optimal model

# include all explanatory variables and as many interactions
# Example: Phenotype ~ X1 + X2 + X1:X2

## Step 2: Find optimal random structure

# Use REML to compare these (nested) models (with ML the variance estimates are biased).
mod1  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (1 | Individual),  data = df, REML = TRUE)
mod2  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (1 | Individual) + (-1 + X1 | Individual), data = df, REML = TRUE)
mod3  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (X1 | Individual), data = df, REML = TRUE)

lattice::dotplot(lme4::ranef(mod3, condVar = TRUE, whichel = "Individual"))

model_selection(mod1, mod2, mod3)

model_selection <- function(...){
	
	return(data.frame("df"    = AIC(...)$df,
										"AIC"   = AIC(...)$AIC,
										"wAIC"  = round(MuMIn::Weights(AIC(...)),3),
										"AICc"  = MuMIn::AICc(...)$AICc,
										"wAICc" = round(MuMIn::Weights(MuMIn::AICc(...)),3)))
	
}

## Step 3: Find optimal fixed structure

# Use ML and not REML
mod4  <- lme4::lmer(Phenotype ~ X1 + (X1 | Individual), data = df, REML = FALSE)
mod5  <- lme4::lmer(Phenotype ~ X1 + X2 + (X1 | Individual), data = df, REML = FALSE)
mod6  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (X1 | Individual), data = df, REML = FALSE)
	
model_selection(mod4, mod5, mod6)

## Step 4: Present the final model using REML estimation
modf  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (X1 | Individual), data = df, REML = TRUE)
summary(modf)


#####################
# RESIDUAL ANALYSIS #
#####################

# Loy, A., & Hofmann, H. (2013). Diagnostic tools for hierarchical linear models. 
# Wiley Interdisciplinary Reviews: Computational Statistics, 5(1), 48–61.

# "HLMdiag" Package #
# Loy, A., & Hofmann, H. (2014). HLMdiag: A Suite of Diagnostics for Hierarchical Linear Models in R Adam. 
# Journal Of Statistical Software, 56(5), 1–28.

# Upward residual analysis

### Level 1 (conditional) residuals

resid1 <- HLMdiag::HLMresid(modf, level = 1, type = "LS", standardize = "semi")
resid1 <- cbind(df$Individual, resid1)
names(resid1)[1] <- "Individual"
head(resid1)

## check for linearity
qplot(x =  fitted, y = LS.resid, data = resid1, geom = c("point", "smooth")) + ylab("LS level-1 residuals")
qplot(x =  X1,     y = LS.resid, data = resid1, geom = c("point", "smooth")) + ylab("LS level-1 residuals")
qplot(x =  X2,     y = LS.resid, data = resid1, geom = c("point", "smooth")) + ylab("LS level-1 residuals")

## check for homoscedasticity
qplot(x =  fitted, y = semi.std.resid, data = resid1, geom = c("point", "smooth")) + ylab("semi-standardized residuals")
qplot(x =  X1,     y = semi.std.resid, data = resid1, geom = c("point", "smooth")) + ylab("semi-standardized residuals")
qplot(x =  X2,     y = semi.std.resid, data = resid1, geom = c("point", "smooth")) + ylab("semi-standardized residuals")

# Roughly constant interquartile ranges between groups
boxplot(LS.resid ~ Individual, data = resid1)

## Normality (quantile plot of the sem-=standardized level-1 residuals)
HLMdiag::ggplot_qqnorm(x = resid1$semi.std.resid, line = "rlm")

### Level 2 (random effects) residuals

resid2 <- HLMdiag::HLMresid(object = modf, level = "Individual")
head(resid2)

## Normality of the random intercept
hist(resid2$'(Intercept)')
HLMdiag::ggplot_qqnorm(x = resid2$'(Intercept)', line = "rlm", main = "Random intercept")

## Normality of the random slope
hist(resid2$X1)
HLMdiag::ggplot_qqnorm(x = resid2$X1, line = "rlm", main = "Random slope")


######################
# INFLUENCE ANALYSIS #
######################

# Start with diagnotics for the variance components as diagnostics 
# for fixed effects require a specific covariance matrix.

### Diagnostics for variance components

# The relative variance change (RVC is close to 0 when deleted element is not influential)
rvc1 <- HLMdiag::rvc(modf) # Level-1 deletion 
rvc2 <- HLMdiag::rvc(modf, group = "Individual") # Level-2 (Individual-level) deletion

head(rvc1)
# sigma2: residual variance
# D11: random intercept variance
# D22: random slope variance
# D12: covariance between random intercept and slope.

# Visualize influential elements
HLMdiag::dotplot_diag(x = rvc1[ , "D11"], cutoff = "internal", name = "rvc", modify = "dotplot") + ylab("Relative random intercept variance change") + xlab("Value")
HLMdiag::dotplot_diag(x = rvc2[ , "D11"], cutoff = "internal", name = "rvc", modify = "dotplot") + ylab("Relative random intercept variance change") + xlab("Individual")


### Diagnostics for fixed effects

## Changes in parameter values

# Cook's distance (larger values indicate higher leveles of influence)
cooksd1 <- cooks.distance(modf) # Level-1 deletion 
cooksd2 <- cooks.distance(modf, group = "Individual") # Level-2 (Individual-level) deletion

# Visualize influential elements
HLMdiag::dotplot_diag(x = cooksd1, cutoff = "internal", name = "cooks.distance", modify = "dotplot") + ylab("Cook's distance") + xlab("Value")
HLMdiag::dotplot_diag(x = cooksd2, cutoff = "internal", name = "cooks.distance", modify = "dotplot") + ylab("Cook's distance") + xlab("Individual")

# Access to the difference in the parameters associated with the deletion of one element
beta_cdd <- as.numeric(attr(cooksd1, "beta_cdd")[[58]])
names(beta_cdd) <- names(lme4::fixef(modf))
beta_cdd

## Precision of fixed parameters

# Covariance ratio (Close to 1 when the deleted element is not influential)
covratio1 <- HLMdiag::covratio(modf) # Level-1 deletion 
covratio2 <- HLMdiag::covratio(modf, group = "Individual") # Level-2 (Individual-level) deletion

# Visualize influential elements
HLMdiag::dotplot_diag(x = covratio1, cutoff = "internal", name = "covratio", modify = "dotplot") + ylab("Covariance ratio") + xlab("Value")
HLMdiag::dotplot_diag(x = covratio2, cutoff = "internal", name = "covratio", modify = "dotplot") + ylab("Covariance ratio") + xlab("Individual")


### Diagnostics for fitted values

# Leverage: can be defined as the rate of change in the fitted response 
# with respect to the observed response
leverage1 <- HLMdiag::leverage(modf, level = 1)
leverage2 <- HLMdiag::leverage(modf, level = "Individual")

head(leverage1)
# overall: averall leverage
# fixef: the fixed effects leverage
# ranef: the random effects leverage (confounded: depends on the leverage associated with the fixed effects)
# ranef.uc: the unconfounded random effects leverage (modified ranef)

# Visualize influential elements
HLMdiag::dotplot_diag(x = leverage1[ , "overall"], cutoff = "internal", name = "leverage", modify = "dotplot") + ylab("Overall leverage") + xlab("Value")
HLMdiag::dotplot_diag(x = leverage2[ , "overall"], cutoff = "internal", name = "leverage", modify = "dotplot") + ylab("Overall leverage") + xlab("Individual")



# ### jackknife per point
# res    <- jackknife_point(modf, df)
# DBfit  <- res$DBfit
# DBbeta <- res$DBbeta
# 
# ggplot(data=DBfit, aes(x=X1, y=fitted, color = as.factor(rank))) +
# 	geom_point() + 
# 	stat_smooth(method = "lm", alpha = 0) + 
# 	theme(legend.position="none")
# 
# 
# ### jackknife per group
# res    <- jackknife.group(mod, df)
# DBfit  <- res$DBfit
# DBbeta <- res$DBbeta
# 
# ggplot(data=DBfit, aes(x=X1, y=fitted, color = as.factor(rank))) +
# 	geom_point() + 
# 	stat_smooth(method = "lm", alpha = 0) + 
# 	theme(legend.position="none")


###################
# AUTOCORRELATION #
###################


acf(residuals(modf))
#variogram()





########################################################################
################### Extract estimated parameters #######################
########################################################################

lme4::fixef(modf)    # fixed effect coefficients
arm::se.fixef(modf) # standard error of fixed effect coefficients

as.data.frame(lme4::VarCorr(modf)) # get random effect (variances)
arm::se.ranef(modf)  # standard error of fixed and random effect coefficients

ran_effect <- lme4::ranef(modf, condVar = TRUE, whichel = "Individual") # random effect for each individual
lattice::dotplot(ran_effect) # plot random effect for each Individual with the standard error

########################################################################
####################### PARAMETRIC BOOTSTRAP ###########################
########################################################################

mySumm <- function(.) {
	c(beta0 = lme4::fixef(.)["(Intercept)"],            # Intercept fixed effect coefficients
		beta1 = lme4::fixef(.)["X1"],                     # Slope fixed effect coefficients
		sig01 = as.data.frame(lme4::VarCorr(.))$vcov[1],  # Intercept random effect variance
		sig02 = as.data.frame(lme4::VarCorr(.))$vcov[1],  # Slope random effect variance
		sigma = sigma(.)^2                          # residual variance
	)
}

boo01 <- lme4::bootMer(modf, mySumm, nsim = 100, .progress = "txt", seed = 101)

## intercept (fixed effect)
lme4::fixef(modf)["(Intercept)"]
(beta0 <- boot::boot.ci(boo01, index=1, type=c("norm", "basic", "perc")))# beta0

## slope (fixed effect)
lme4::fixef(modf)["X1"]
(beta1 <- boot::boot.ci(boo01, index=2, type=c("norm", "basic", "perc")))# beta1

## intercept standard deviation (random effect)
as.data.frame(lme4::VarCorr(modf))$vcov[1]
(sig01 <- boot::boot.ci(boo01, index=3, type=c("norm", "basic", "perc")))# sig01

## slope standard deviation (random effect)
as.data.frame(lme4::VarCorr(modf))$vcov[4]
(sig02 <- boot::boot.ci(boo01, index=4, type=c("norm", "basic", "perc")))# sig02

## residual standard deviation
sigma(mod)^2
(sigma <- boot::boot.ci(boo01, index=5, type=c("norm", "basic", "perc")))# sigma

#########################################
