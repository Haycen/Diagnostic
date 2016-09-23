# load packages and functions
source("diagnostic_fcns.r")
library(lme4)
library(car)
library(lattice)
library(arm)
library(ggplot2)

## Load dataset
df <- read.csv("datasets/dataset_1.csv")

########################################################################
############ Issues to consider before fitting a model #################
########################################################################

#### Data distribution
distribution.plot(df$Phenotype)


#################################################################################################
#### Check for missing data
df$missing <- ifelse(is.na(df$Phenotype), TRUE, FALSE)

# test if predictors of missing data are different from the not missing data
t.test(X1~missing, df)
boxplot(X1~missing, df)


#################################################################################################
#### OUTLIERS

# Response variable
par(mfrow = c(2, 2))
# Response variable
boxplot(df$Phenotype, ylab = "Phenotype", main = "Boxplot")
dotchart(df$Phenotype,
				 ylab = "Order of observations",
				 xlab = "Phenotype", main = "Cleveland dotplot")
boxplot(df$Phenotype~df$Individual, ylab = "X1", xlab = "Individual", main = "Boxplot")
dotchart(df$Phenotype,
				 groups = factor(df$Individual),
				 ylab = "Individual", xlab = "Phenotype",
				 main = "Cleveland dotplot", pch = df$Individual)


# Explanatory variable
par(mfrow = c(2, 2))
# Explanatory variable
boxplot(df$X1, ylab = "X1", main = "Boxplot")
dotchart(df$X1, ylab = "Order of observations", xlab = "X1", main = "Cleveland dotplot")
boxplot(df$X1~df$Individual, ylab = "X1", xlab = "Individual", main = "Boxplot")
dotchart(df$X1,
				 groups = factor(df$Individual),
				 ylab = "Individual", xlab = "X1",
				 main = "Cleveland dotplot", pch = df$Individual)



par(mfrow=c(1,1))
mod_lm  <- lm(Phenotype ~ X1, data = df)

car::leveragePlots(mod_lm) # leverage plots
car::outlierTest(mod_lm) # Bonferonni p-value for most extreme obs
car::influencePlot(mod_lm)
plot(mod_lm) # display diagnostic plot for an lm object


# Leverage tests (Laszlo)
lev = hat(model.matrix(mod_lm))
max(stats:::influence(mod_lm)$hat)
level1=2*(length(coefficients(mod_lm))+1)/length(residuals(mod_lm)) #laverage should be smaller than this
level2=3*(length(coefficients(mod_lm))+1)/length(residuals(mod_lm)) #or this

par(mfrow=c(1,2))
plot(lev, ylim=c(0, level2))
abline(h=level1, col="red", lty=2, lwd=1)
abline(h=level2, col="red", lty=1, lwd=2)
cutoff <- 4/((nrow(df)-length(mod_lm$coefficients)-2))
plot(mod_lm, which=4, cook.levels=cutoff)


#### Influencial points (using influenceME package)
library(influence.ME)

inl.diag.obs <- influence.ME:::influence(mod, obs = T)

#cooks.distance
inl.diag.cooks <- cooks.distance(inl.diag.obs)
plot(inl.diag.obs, which="cook",
		 cutoff=.17, sort=TRUE,
		 xlab="Cook´s Distance",
		 ylab="Individual's ID")

# dfbetas
inl.diag.dfB <- dfbetas.estex(inl.diag.obs, sort = TRUE, to.sort = "X1")
plot(inl.diag.obs, which="dfbetas",
		 cutoff=.17,sort = TRUE, to.sort = "X1",
		 xlab="DFBETAS",
		 ylab="Individuals ID")
plot(inl.diag.obs, which="dfbetas",
		 cutoff=.17,sort = TRUE, to.sort = "(Intercept)",
		 xlab="DFBETAS",
		 ylab="Individuals ID")

### jackknife per point
res    <- jackknife_point(mod, df)
DBfit  <- res$DBfit
DBbeta <- res$DBbeta

ggplot(data=DBfit, aes(x=X1, y=fitted, color = as.factor(rank))) +
	geom_point() + 
	stat_smooth(method = "lm", alpha = 0) + 
	theme(legend.position="none")

#### Influencial groups (using influenceME package)

# visualize the data per individual
plot(df$Individual, df$Phenotype)

inl.diag.group <- influence.ME:::influence(mod, group = "Individual")

#cooks.distance
inl.diag.cooks <- cooks.distance(inl.diag.group)
plot(inl.diag.group, which="cook",
		 cutoff=.17, sort=TRUE,
		 xlab="Cook´s Distance",
		 ylab="Individual's ID")

inl.diag.dfB <- dfbetas.estex(inl.diag.group, sort = TRUE, to.sort = "X1")
plot(inl.diag.group, which="dfbetas",
		 cutoff=.17,sort = TRUE, to.sort = "X1",
		 xlab="DFBETAS",
		 ylab="Individuals ID")
plot(inl.diag.group, which="dfbetas",
		 cutoff=.17,sort = TRUE, to.sort = "(Intercept)",
		 xlab="DFBETAS",
		 ylab="Individuals ID")

### jackknife per group
res    <- jackknife.group(mod, df)
DBfit  <- res$DBfit
DBbeta <- res$DBbeta

ggplot(data=DBfit, aes(x=X1, y=fitted, color = as.factor(rank))) +
	geom_point() + 
	stat_smooth(method = "lm", alpha = 0) + 
	theme(legend.position="none")


#################################################################################################

















#### Collinearity between the model predictors

# Correlation between two predictors
with(df,cor(X1, X2))

# VIF
df   <- df[!is.na(df$Phenotype),] # Remore rows with missing values
mod  <- lmer(Phenotype ~ 1 + X1 + X2 + X1X2 + (X1|Individual), data = df)
vif.mer(mod)

########################################################################
########################### Fit the model ##############################
########################################################################
df   <- df[!is.na(df$Phenotype),] # Remore rows with missing values

# fit mixed-effect model
mod  <- lmer(Phenotype ~ 1 + X1 + (X1|Individual), data = df)

summary(mod)

#### examining residuals #### 
# centered around 0
# normality 
# homogeneity of variance
diagnostics.plot(mod, df)

# Check for residual pattern within individuals and difference between individuals      
lattice::xyplot(residuals(mod) ~ fitted(mod) | df$Individual,
			 main = "Residual pattern within Individuals",
			 panel=function(x, y){ 
			 	panel.xyplot(x, y) 
			 	panel.loess(x, y, span = 0.75) 
			 	panel.lmline(x, y, lty = 2)  # Least squares broken line
			 } 
)

### Check for autocorrelation in the residuals (independency)
acf(residuals(mod))

#### Test for outliers
mod_lm  <- lm(Phenotype ~ 1 + X1, data = df)
outlierTest(mod_lm) # Bonferonni p-value for most extreme obs
leveragePlots(mod_lm) # leverage plots
plot(mod_lm) # display diagnostic plot for an lm object



########################################################################
################### Extract estimated parameters #######################
########################################################################

fixef(mod)    # fixed effect coefficients
se.fixef(mod) # standard error of fixed effect coefficients

as.data.frame(VarCorr(mod)) # get random effect (variances)
se.coef(mod)  # standard error of fixed and random effect coefficients

ran_effect <- ranef(mod, condVar=TRUE, whichel = "Individual") # random effect for each individual
lattice::dotplot(ran_effect) # plot random effect for each Individual with the standard error

########################################################################
####################### PARAMETRIC BOOTSTRAP ###########################
########################################################################
require(boot)

mySumm <- function(.) {
	c(beta0 = fixef(.)["(Intercept)"],            # Intercept fixed effect coefficients
		beta1 = fixef(.)["X1"],                     # Slope fixed effect coefficients
		sig01 = as.data.frame(VarCorr(.))$vcov[1],  # Intercept random effect variance
		sig02 = as.data.frame(VarCorr(.))$vcov[1],  # Slope random effect variance
		sigma = sigma(.)^2                          # residual variance
	)
}

boo01 <- bootMer(mod, mySumm, nsim = 100, .progress = "txt", seed=101)

## intercept (fixed effect)
fixef(mod)["(Intercept)"]
(beta0 <- boot.ci(boo01, index=1, type=c("norm", "basic", "perc")))# beta0

## slope (fixed effect)
fixef(mod)["X1"]
(beta1 <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc")))# beta1

## intercept standard deviation (random effect)
as.data.frame(VarCorr(.))$vcov[1]
(sig01 <- boot.ci(boo01, index=3, type=c("norm", "basic", "perc")))# sig01

## slope standard deviation (random effect)
as.data.frame(VarCorr(.))$vcov[4]
(sig02 <- boot.ci(boo01, index=4, type=c("norm", "basic", "perc")))# sig02

## residual standard deviation
sigma(mod)^2
(sigma <- boot.ci(boo01, index=5, type=c("norm", "basic", "perc")))# sigma

#########################################
