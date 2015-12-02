diagnostics.plot<-function(mod.res, data){
  
	old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 1, 0.5))
  
  # Histogram of residuals
  hist(residuals(mod.res), probability=TRUE, xlab="", ylab="", main="")
  abline(v=0, lwd=2, lty=2, col = "red")
  mtext(text="histogram of residuals", side=3, line=0)
  
  # qq-plot of residuals
  x <- seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(residuals(mod.res))))
  qqnorm(residuals(mod.res), main="", pch=19)
  qqline(residuals(mod.res), lwd=2, lty=2)
  mtext(text="qq-plot of residuals", side=3, line=0)
  
  # residuals against fitted values
  plot(fitted(mod.res), residuals(mod.res), pch=19)
  lines(smooth.spline(fitted(mod), residuals(mod)), col = "blue",lwd=2)
  abline(h=0, lty=2)
  mtext(text="residuals against fitted values", side=3, line=0)
  
  # residuals for  each individual
  boxplot(residuals(mod.res) ~ Individual, data=data)
  abline(h=0, lwd=2, lty=2)
  mtext(text="residuals for  each individual", side=3, line=0)
  
  par(old.par)
}

distribution.plot <- function(var){
	par(mfrow=c(2, 2))
	par(mar=c(3, 3, 1, 0.5))
	
	hist(var, probability=TRUE, xlab="", ylab="", main="")
	abline(v=0, lwd=2, lty=2, col = "red")
	mtext(text="histogram of sampled data", side=3, line=0)
	
	qqnorm(var, main="", pch=19)
	qqline(var, lwd=2, lty=2)
	mtext(text="qq-plot of sampled data", side=3, line=0)
	
	#how random distribution should look
	r.test <- rnorm(n=length(var), mean=mean(var, na.rm = T), sd=sd(var, na.rm = T)) 
	
	hist(r.test, probability=TRUE, xlab="", ylab="", main="")
	abline(v=0, lwd=2, lty=2, col = "red")
	mtext(text="histogram of theorical data", side=3, line=0)
	
	qqnorm(r.test, main="", pch=19)
	qqline(r.test, lwd=2, lty=2)
	mtext(text="qq-plot of theorical data", side=3, line=0)
}

lev.thresh<-function(model.res){
	k=length(coefficients(model.res))
	n=length(residuals(model.res))
 return(2*(k+1)/n)
}

overdisp.test<-function(x){
  pr=residuals(x, type ="pearson")
  sum.dp=sum(pr^2)
  if(class(x)[[1]]=="mer"){
    if(length(grep(x=x@call, pattern="poisson"))==0){
      stop("oops, this isn't a model with poisson family... this function doesn't make sense")
    }
    xdf=length(residuals(x))-length(fixef(x))
  }else{
    if (x$family[[1]]!="poisson"){
      stop("oops, this isn't a model with poisson family... this function doesn't make sense")
    }
    xdf=length(residuals(x))-length(x$coefficients)
  }
  return(data.frame(chisq=sum.dp, df=xdf, P=1-pchisq(sum.dp, xdf), dispersion.parameter=sum.dp/xdf))
}

vif.mer <- function (fit) {
	## adapted from rms::vif
	
	v <- vcov(fit)
	nam <- names(fixef(fit))
	
	## exclude intercepts
	ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
	if (ns > 0) {
		v <- v[-(1:ns), -(1:ns), drop = FALSE]
		nam <- nam[-(1:ns)]
	}
	
	d <- diag(v)^0.5
	v <- diag(solve(v/(d %o% d)))
	names(v) <- nam
	v
}

jackknife_point <- function(mod, df){
	
	DBfit  <- c()
	DBbeta <- c()
	
	for (i in 1:nrow(df)){
		red    <- update(mod, data=df[-i,])
		DBfit  <- rbind(DBfit, data.frame("fitted" = fitted(red), 
																			"X1" = df$X1[-i],
																			"rank" = i))
		DBbeta <- rbind(DBbeta, data.frame("intercept" = coef(red)$Individual[,1], 
																			 "slope"     = coef(red)$Individual[,2],
																			 "rank" = i))
	}
	
	return(list("DBfit" = DBfit, "DBbeta" = DBbeta))
}

jackknife.group <- function(mod, df){
	
	DBfit  <- c()
	DBbeta <- c()
	
	for (i in unique(df$Individual)){
		
		new_data <- df[df$Individual != i,]
		red    <- update(mod, data=new_data)
		DBfit  <- rbind(DBfit, data.frame("fitted" = fitted(red), 
																			"X1" = new_data$X1,
																			"rank" = i))
		DBbeta <- rbind(DBbeta, data.frame("intercept" = coef(red)$Individual[,1], 
																			 "slope"     = coef(red)$Individual[,2],
																			 "rank" = i))
	}
	
	return(list("DBfit" = DBfit, "DBbeta" = DBbeta))
}

