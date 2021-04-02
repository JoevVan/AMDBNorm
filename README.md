# AMDBNorm
This is a new algorithm, adjustment mean DBNorm (AMDBNorm), which is based on a probability distribution to correct batch effects while preserving biological variation. We expect that AMDBNorm will improve statistical analysis ability for some rare diseases and help with disease prediction, new subtype discovery and biomarker researches.
# Install in R
install.packages("devtools")

library("devtools")

install_github("JoevVan/AMDBNor")

PS: if lazy loading failed to install a required package "DBNorm", "future.apply", "future", "reshape2", please install it manually before running AMDBNorm.

devtools::install_github("mengqinxue/DBNorm")

install.package("future.apply")

install.package("future")

install.package("reshape2")

library('DBNorm')

library('reshape2')

install.package("future.apply")

# Loading build-in datasets
#These test datasets are part of the real datasets, including three batches of samples(‘Batch_1’, ‘Batch_2’, ‘Batch_3’, ‘Group_1’, ‘Group_2’, ‘Group_3’).

data(Test_data)
# Define reference distribution
Batch_1[1:4,1:4]#Reference Batch

dis <- genDistData(melt(Batch_1)[,3], 500)
# Visualising distribution datasets.
visDistData(dis, "P", "Reference distribution", "Range", "Probability")
# Using POLT function to fit reference distribution.
fit <- polyFit(dis,9)# Visualising fitting results

visFitting(fit, "Reference distribution", "Range", "Probability")
# Batch correction using AMDBNorm.
AMDB_Batch_2 <- AMDBNorm(data = Batch_2,ref_batch = Batch_1,fit = fit)

AMDB_Batch_1 <- AMDBNorm(data = Batch_1,ref_batch = Batch_1,fit = fit)

#AMDB_Batch_2 <- AMDBNorm(data = Batch_2,ref_Batch = Batch_1,fit = fit,thread = TRUE)

#AMDB_Batch_1 <- AMDBNorm(data = Batch_1,ref_Batch = Batch_1,fit = fit,thread = TRUE)
# Application of AMDBNorm in unbalanced experimental design.
#Suppose the experimental design has two biological conditions: normal(N) and tumor(T)
# Define reference distributions
dis_N <- genDistData(melt(Batch_1[,Group_1=='normal'])[,3], 500)

dis_T <- genDistData(melt(Batch_1[,Group_1=='tumor'])[,3], 500)
# visualising reference distributions
visDistData(dis_N, "P", "Reference distribution (N)", "Range", "Probability")

visDistData(dis_T, "P", "Reference distribution (T)", "Range", "Probability")
# Using POLT function to fit reference distributions
fit_N <- polyFit(dis_N,9)

fit_T <- polyFit(dis_T,9)
# Batch correction using AMDBNorm.
AMDB_Batch_3N <- AMDBNorm(data = Batch_3[,Group_3=='normal'],
                          ref_batch = Batch_1[,Group_1=='normal'], fit = fit_N)
                          
AMDB_Batch_3T <- AMDBNorm(data = Batch_3[,Group_3=='tumor'],
                          ref_batch = Batch_1[,Group_1=='tumor'], fit = fit_T)
# The result of Batch correction of Batch 3.
AMDB_Batch_3 <- cbind(AMDB_Batch_3N,AMDB_Batch_3T)

Group_A3 <- c(rep('normal',ncol(AMDB_Batch_3N)),
              rep('tumor',ncol(AMDB_Batch_3T)))
# Batch correction result
AMDB_Batch1_2_3 <- cbind(AMDB_Batch_1,AMDB_Batch_2,AMDB_Batch_3)

AMDB_Batch1_2_3[1:4,1:4]

boxplot(AMDB_Batch1_2_3)

Group <- c(Group_1,Group_2,Group_A3)

whichBatch <- c(rep('Batch_1',ncol(Batch_1)),rep('Batch_2',ncol(Batch_2)),
                rep('Batch_3',ncol(Batch_3)))
