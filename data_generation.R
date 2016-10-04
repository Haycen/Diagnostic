library(data.table)

#############
# Dataset 1 # (used as example)
#############

# Model: Phenotype ~ X1 + X2 + X1:X2 + (1 | Individual) + (-1 + X1 | Individual)

# Read data
orig_data_1 <- fread("./datasets/dataset_original_1.csv")
# Subset data
dataset_1 <- orig_data_1[ , .(Individual, Time, Phenotype, X1, X2)]

# Parameters:
# NI      = 12
# NR      = 15
# d1_b0   = 0.5
# d1_b1   = 0.8
# d1_b2   = 0.2
# d1_b12  = 0.5
# d1_VI   = 1.0
# d1_VS1  = 0.5
# d1_Ve   = 0.05

# Add missing data (15% of response variable)
dataset_1[sample(nrow(dataset_1), round(nrow(dataset_1) * 0.15), replace = FALSE), Phenotype := NA]

write.csv(x = dataset_1, file = "./datasets/dataset_1.csv", row.names = FALSE)


#############
# Dataset 2 # Collinearity + individual outlier
#############

# Model: Phenotype ~ X1 + X2 + (1 | Individual) + (-1 + X1 | Individual) + (-1 + X2 | Individual)
# X3 is collinear to X2
# Individual 21 is an outlier

# Read data
orig_data_2 <- fread("./datasets/dataset_original_2.csv")
# Subset data
dataset_2 <- orig_data_2[ , .(Individual, Time, Phenotype, X1, X2)]

ind_21 <- dataset_2[Individual == 14, ]
ind_21[ , ':='(Individual = 21, Phenotype = Phenotype + rnorm(.N, 5, 0.3))]

dataset_2 <- rbind(dataset_2, ind_21)
dataset_2[ , X3 := 0.7 * X2 + rnorm(.N, 0, 0.3)]

# Parameters:
# NI      = 20 
# NR      ~ 10
# d1_b0   = 0.9
# d1_b1   = 1.2
# d1_b2   = 0.2
# d1_VI   = 1.0
# d1_VS1  = 0.6
# d1_VS2  = 0.2
# d1_Ve   = 0.05

write.csv(x = dataset_2, file = "./datasets/dataset_2.csv", row.names = FALSE)


#############
# Dataset 3 # X1 has a left-skewed distribution (need transformation)
#############

# Model: Phenotype ~ log(X1) + (1 | Individual)


# Read data
orig_data_3 <- fread("./datasets/dataset_original_3.csv")
orig_data_3[ , X1 :=  exp(X1)]

# Subset data
dataset_3 <- orig_data_3[ , .(Individual, Time, Phenotype, X1)]

# rgamma(1000, 3)

# Parameters:
# NI      = 25 
# NR      = 10
# d1_b0   = 3.4
# d1_b1   = 1.3
# d1_VI   = 0.9
# d1_Ve   = 0.05

write.csv(x = dataset_3, file = "./datasets/dataset_3.csv", row.names = FALSE)


#############
# Dataset 4 # 
#############

# Model: Phenotype ~ X1 + X2 + X1:X2 + (X1 | Individual) + (-1 + X1:X2 | Individual)

# Read data
orig_data_4 <- fread("./datasets/dataset_original_4.csv")
# orig_data_4[ , Phenotype :=  0 + I + (B1 * X1) + rgamma(.N, 1)]

# Subset data
dataset_4 <- orig_data_4[ , .(Individual, Time, Phenotype, X1, X2)]

# Parameters:
# NI        = 25 
# NR        ~ 10
# b0        = 3.4
# b1        = 1.3
# b2        = 0
# b12       = 3
# VI        = 0.9
# VS1       = 0.9
# Corr I-S1 = 0.8
# VS2       = 0
# VS12      = 0.5
# Ve        = 0.05

write.csv(x = dataset_4, file = "./datasets/dataset_4.csv", row.names = FALSE)


#############
# Dataset 5 # Residual variance increase qith X1
#############

# Model: Phenotype ~ X1 + X2 + (1 | Individual)

# Read data
orig_data_5 <- fread("./datasets/dataset_original_5.csv")

orig_data_5[ , X_tmp := X1 - min(X1)]
orig_data_5[ , me := rnorm(.N,0,X_tmp*2)]
orig_data_5[ , Phenotype := B0 + I + (B1 * X1) + (B2 * X2) + me]

# Subset data
dataset_5 <- orig_data_5[ , .(Individual, Time, Phenotype, X1, X2)]

# Parameters:
# NI   = 10 
# NR   ~ 9
# X1   = cyclic
# b0   = 0
# b1   = 0.5
# b2   = 0.3
# VI   = 0.7
# Ve   = 0.05

write.csv(x = dataset_5, file = "./datasets/dataset_5.csv", row.names = FALSE)
