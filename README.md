# Outline

The document contains five R code documents and one pdf showing a simplified simulation
result of original paper of P15. The following describes functionalities of different documents, and give instructions on how to use them. 

Be aware: QUIC, mvtnorm, caret, pracma, parallel packages are required. Some packages might require users have R version above 3.5.3. In the first line of Main.R, users have to reset the default working directory. 

User manual: Users do not have to do any modification in Sampling.R, Estimation.R and Precision.R. But in Crossvalidation.R, the users must reset the number of cores to allow parallel computing in line 16 according to how many cores his/her computers have, otherwise the program is doomed to crash down. To know how many cores you have, try "require(parallel)" and "detectCores()". In Main.R, users can change experimental variables in line 6, where n is sample size, B is number of replications (strongly recommended to set to 1 if you want a quick result), p is number of dimension, mode=1-4 corresponds to four different types precision matrices sampling schemes. In line 28, for(m in 1:4) can be changed to 1:1 if users do not have enough time to run the codes. Performance measurements will be given in list "dat" when the program terminates.

Computation speed reference: Under setting n=200, p=120, B=1, mode=1, m=1,2 needs 4mins, m=3,4 needs 40mins. 

## 1. Sampling.R

This document contains a function (Sampling) that attempts to realize the sampling scheme proposed in P15 paper. This function consists of two arguments, where n indicates the number of total sample points, defaulted by n=200 as in P15, Omega is the type of precision matrix used to generate "signal" part of proposed model. The function will return a list of seven different types of data sets according to different contamination schemes in the paper, including clean, rowwise, cellwise contamination of Gaussian models, as well as heavier-tail models.

## 2. Precision.R

The paper considers four different types of precision matrix. This document contains a function (Precision) that aims to automatically generate desired matrices. Argument p is dimension of matrix, defaulted by 120; mode is equal to 1-4, where 1 generates Banded type, 2 Sparse type, 3 Dense type and 4 Diagonal type. The output of this function can be directly fed into Sampling function in 1.

## 3. Estimation.R

The document contains a function (Estimation) that returns estimated covariance & precision matrix given data X. When Inv=T, given estimated covariance according to different specified methods, it returns estimated Precision matrix based on Graphical Lasso. 

Argument method consists of "SampleCov", "SpearmanU", "Spearman" and "Kendall", "InvCov", their meaning is the same as in original paper; lambda is the penalty parameter; Inv is a logical value indicating whether it is necessary to return estimated Precision matrix.

## 4. Crossvalidation.R 

This document contains function (Crossvalidation) that returns the best value of lambda used to penalize Graphical Lasso. The choice of values of lambda is described in P15 paper.  Argument fold indicates how many folds you want to use in crossvalidation selection; method has the same meaning as in 3.

## 5. Main.R

This document consists of codes used to conduct simulation study.










