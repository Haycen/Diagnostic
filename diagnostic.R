# load packages and functions
# source("diagnostic_fcns.r")
library(lme4)
library(car)
library(lattice)
library(arm)
library(GGally)
library(ggplot2)
library(HLMdiag)
library(boot)
library(BaylorEdPsych)
library(psych)
library(MuMIn)
# install.packages("mvnmle")

## Load dataset
df <- read.csv("datasets/dataset_2.csv")

head(df)
str(df)

df$Individual <- as.factor(df$Individual)

# Examine sampling design
table(df$Individual)

# Visualize sampling time per individual
ggplot(data = df, aes(x = Time,
											y = Individual,
											color = Individual)) +
	geom_point() +
	ggtitle("Sampling time per individual") +
	xlab("Time") +
	ylab("Individuals")


########################################################################
############ Issues to consider before fitting a model #################
########################################################################

#### Data overview

GGally::ggpairs(
	df, columns = c("X1", "X2", "Phenotype"),
	diag  = list(continuous = "barDiag", colour = "grey"),
	lower = list(
		continuous = "smooth",
		mapping    = aes(color = Individual)
	)
)


################
# MISSING DATA #
################

# Nakagawa, S., & Freckleton, R. P. (2008). Missing inaction: the dangers of ignoring missing data. 
# Trends in Ecology and Evolution, 23(11), 592–596.

# Visualize missing data
Amelia::missmap(df[ , c("Phenotype", "X1", "X2")])

# missing completely at random (MCAR): the probability of missing
# data in one variable is not related to any other variable in the data set

# missing at random (MAR): the probability of missing data in a variable
# is related to some other variable(s) in the data set

# missing not at random (MNAR):the probability of missing data in a variable
# is associated with this variable itself, even after controlling 
# for other observed (related) variables

## MCAR vs MAR and/or MNAR

# test if predictors of missing data are different from the not missing data
df$missing <- ifelse(is.na(df$Phenotype), 1, 0)
t.test(X1~missing, df)
boxplot(X1~missing, df)

t.test(X2~missing, df)
boxplot(X2~missing, df)

# Little, R. J. A. (1988). A test of missing completely at random for multivariate data 
# with missing values. Journal of the American Statistical Association, 83(404), 1198–1202.

# Little's MCAR test
res <- BaylorEdPsych::LittleMCAR(df[ , c("Phenotype", "X1", "X2")])
# If significant data set contains non-MCAR missingness
res$chi.square
res$p.value


# create the missingness matrix
Missingness <- ifelse(is.na(df[ , c("Phenotype", "X1", "X2")]) == TRUE, 1, 0)

# combine the original dataset with the missingness matrix
MissData <- data.frame(df[ , c("Phenotype", "X1", "X2")], Missingness)
psych::pairs.panels(MissData, ellipses = FALSE, method = "spearman")


## Methods for Missing data
# Data deletion
# Single imputation
# Multiple imputation
# Data augmentation

# Book:
# Fox, G.A., Negrete-Yankelevich, S., Sosa, V.J. (2015). Ecological Statistics: 
# Contemporary Theory and Application. Oxford University Press.



################
# COLLINEARITY #
################

# VIF: Variance inflation factor
# if VIF < threshold, no collinearity
# Threshold =  10 (Montgomery, D.C. & Peck, E.A. 1992. Wiley)
# Threshold =  3 (Zuur, A.F. et al. 2010. Methods in Ecology and Evolution)

mod  <- lm(Phenotype ~ X1 + X2, data = df)
car::vif(mod)

# Collinearity can be solved by dropping collinear covariates:
# Using VIF or (prehaps better) use common sens and biological knowledge



###################
# MODEL SELECTION #
###################

# Zuur, A. F. et al. (2009). Mixed Effects Models and Extensions in Ecology with R.

### Step 1: beyond optimal model

# include all explanatory variables and as many interactions
# Example: Phenotype ~ X1 + X2 + X1:X2

### Step 2: Find optimal random structure

# Use REML to compare these (nested) models (with ML the variance estimates are biased).
mod1  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (1 | Individual),  data = df, REML = TRUE)
mod2  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (1 | Individual) + (-1 + X1 | Individual), data = df, REML = TRUE)
mod3  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (X1 | Individual), data = df, REML = TRUE)


lattice::dotplot(lme4::ranef(mod3, condVar = TRUE, whichel = "Individual"))


data.frame("df"    = AIC(mod1, mod2, mod3)$df,
	 				 "AIC"   = AIC(mod1, mod2, mod3)$AIC,
	 				 "wAIC"  = round(MuMIn::Weights(AIC(mod1, mod2, mod3)),3),
	 				 "AICc"  = MuMIn::AICc(mod1, mod2, mod3)$AICc,
	 				 "wAICc" = round(MuMIn::Weights(MuMIn::AICc(mod1, mod2, mod3)),3))


### Step 3: Find optimal fixed structure

# Use ML and not REML
mod4  <- lme4::lmer(Phenotype ~ X1 + (1 | Individual), data = df, REML = FALSE)
mod5  <- lme4::lmer(Phenotype ~ X1 + X2 + (1 | Individual), data = df, REML = FALSE)
mod6  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (1 | Individual), data = df, REML = FALSE)

data.frame("df"    = AIC(mod4, mod5, mod6)$df,
					 "AIC"   = AIC(mod4, mod5, mod6)$AIC,
					 "wAIC"  = round(MuMIn::Weights(AIC(mod4, mod5, mod6)),3),
					 "AICc"  = MuMIn::AICc(mod4, mod5, mod6)$AICc,
					 "wAICc" = round(MuMIn::Weights(MuMIn::AICc(mod4, mod5, mod6)),3))


### Step 4: Present the final model using REML estimation
modf  <- lme4::lmer(Phenotype ~ X1 + X2 + X1:X2 + (1 | Individual) + (-1 + X1 | Individual), data = df, REML = TRUE)

summary(modf)



#####################
# RESIDUAL ANALYSIS #
#####################

# Response variable
par(mfrow = c(2, 2))
# Response variable
boxplot(df$Phenotype, ylab = "Phenotype", main = "Boxplot")
dotchart(df$Phenotype,
				 ylab = "Order of observations",
				 xlab = "Phenotype", main = "Cleveland dotplot")

boxplot(df$Phenotype~df$Individual, ylab = "Phenotype", xlab = "Individual", main = "Boxplot")
dotchart(df$Phenotype,
				 groups = df$Individual,
				 ylab = "Individual", xlab = "Phenotype",
				 main = "Cleveland dotplot", pch = df$Individual)


# Explanatory variable
par(mfrow = c(2, 2))
# Explanatory variable
boxplot(df$X1, ylab = "X1", main = "Boxplot")
dotchart(df$X1, ylab = "Order of observations", xlab = "X1", main = "Cleveland dotplot")
boxplot(df$X1~df$Individual, ylab = "X1", xlab = "Individual", main = "Boxplot")
dotchart(df$X1,
				 groups = df$Individual,
				 ylab = "Individual", xlab = "X1",
				 main = "Cleveland dotplot", pch = df$Individual)

par(mfrow = c(1, 1))



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

diagnostics.plot(modf, df)


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
HLMdiag::dotplot_diag(x = rvc2[ , "D11"], cutoff = "internal", name = "rvc") + ylab("Relative random intercept variance change") + xlab("Individual")


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
# overall: overall leverage
# fixef: the fixed effects leverage
# ranef: the random effects leverage (confounded: depends on the leverage associated with the fixed effects)
# ranef.uc: the unconfounded random effects leverage (modified ranef)

# Visualize influential elements
HLMdiag::dotplot_diag(x = leverage1[ , "overall"], cutoff = "internal", name = "leverage", modify = "dotplot") + ylab("Overall leverage") + xlab("Value")
HLMdiag::dotplot_diag(x = leverage2[ , "overall"], cutoff = "internal", name = "leverage", modify = "dotplot") + ylab("Overall leverage") + xlab("Individual")


###################
# AUTOCORRELATION #
###################


acf(residuals(modf))

#variogram()



################################
# Extract estimated parameters #
################################

lme4::fixef(modf)   # fixed effect coefficients
arm::se.fixef(modf) # standard error of fixed effect coefficients

as.data.frame(lme4::VarCorr(modf)) # get random effect (variances)
arm::se.ranef(modf)  # standard error of fixed and random effect coefficients

ran_effect <- lme4::ranef(modf, condVar = TRUE, whichel = "Individual") # random effect for each individual
lattice::dotplot(ran_effect) # plot random effect for each Individual with the standard error



########################
# PARAMETRIC BOOTSTRAP #
########################

mySumm <- function(.) {
	c(beta  = lme4::fixef(.),                       # fixed effect coefficients
		ranv  = as.data.frame(lme4::VarCorr(.))$vcov  # random effect variances
	)
}

# Run bootstrap
boo01 <- lme4::bootMer(modf, mySumm, nsim = 100, .progress = "txt", seed = 101)

# Custom function to extract parameter estimate confidence intervals
bCI.tab <- function(b, original, ind=length(b$t0), type="perc", conf=0.95) {
	
	btab0 <- t(sapply(as.list(seq(ind)), function(i) boot::boot.ci(b, index = i, conf = conf, type = type)$percent))
	btab  <- btab0[,4:5]
	rownames(btab) <- names(b$t0)
	a <- (1 - conf)/2
	a <- c(a, 1 - a)
	pct <- stats:::format.perc(a, 3)
	btab <- cbind(btab[,1], original, btab[,2])
	colnames(btab) <- c(pct[1], "original", pct[2])
	
	return(btab)
}

# Calculate confidence intervals
bCI.tab(boo01, mySumm(modf))
