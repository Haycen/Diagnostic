library(squid)

##################################################################
###################### SQuID inputs ##############################
##################################################################

parameters <- data.frame( "NI"        = c(50 , 100), # Number of individuals
                          "NR"        = c(6, 3),   # Number of samples per individual
                          "VI"        = 0.5, # Among-individual variance (intercept)
                          
                          "B0"        = 0.5,   # Population mean phenotype (intercept)
                          "Ve"        = 0.5, # Measurement error variance
                          
                          "NG"        = 1,  # Number of high-level groups
                          "VG"        = 0, # High-level group variance
                          
                          "Tmax"      = 100, # total time steps
                          
                          "ST_ind" = FALSE # If FALSE individuals are not sampled at the same time
                          )

##################################################################
################### Run all simulations ##########################
##################################################################
for (i in 1:nrow(parameters)) { 
  
  inputs <- as.list(parameters[i,])
  inputs$Vind      <- matrix(0, nrow = 4, ncol = 4)
  inputs$Vind[1,1] <- parameters$VI[i]
  inputs$B         <- c(parameters$B0[i], rep(0,3))
  
  dt_tmp <- squid::squidR(inputs)$sampled_data
  dt_tmp <- cbind(dt_tmp, parameters[i,])
  
  if (i == 1) {
    dt <- dt_tmp
  }else {
    dt <- rbind(dt, dt_tmp)
  }
}


# Example of subsetting data 
sub_data <- subset(dt, NI == 100 & NR == 3)
