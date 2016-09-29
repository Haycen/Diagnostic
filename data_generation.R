library(data.table)

#############
# Dataset 1 # (used as example)
#############

# Model: Phenotype ~ X1 + X1 + X2 + X1:X2 (1 | Individual)

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

# Model: Phenotype ~ X1 + X2 + X3 (1 | Individual)
# X3 is collinear to X2

# Read data
orig_data_2 <- fread("./datasets/dataset_original_2.csv")
# Subset data
dataset_2 <- orig_data_2[ , .(Individual, Time, Phenotype, X1, X2)]

dataset_2[ , X3 := 0.7 * X2 + rnorm(.N, 0, 0.3)]
dataset_2[Individual == 16, Phenotype := Phenotype + rnorm(.N, 2, 0.3)]

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









