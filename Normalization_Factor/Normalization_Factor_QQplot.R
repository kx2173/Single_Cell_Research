
library(SPARK)

###############
##__Model 1__##
###############


#### Model 1: Directly generate Row Count with N(mu, sd):


##------------------------------------------------
## Data Generation
##------------------------------------------------
## Generate the countdata based on the info from the realdata and
## patterns summarized from the SpatialDE result

## Fold Change
##-----------------
rm(list = ls())

load(paste0("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_spark.rds")) 

info <- spark@location
tol_counts <- spark@lib_size
# mean(tol_counts) = 29899.95
# sd(tol_counts) = 12990.51

beta <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$coefficients
})

nb <- sapply(1:4, function(x) {
  log(x * exp(median(beta)))
})

tau2 <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$theta[2]
})


#### Model 1 set counts as normal distribution:

set.seed(123) # set seed for reproducibility

nor_count <- round(rnorm(length(tol_counts), mean(tol_counts), 10))


#### Model 1 generate 100 SE and 900 non-SE genes:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35

for (ifc in 3:3) {
  newN <- nor_count 
  for (ipt in 2:2) {
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    numSignal <- 100
    numGene <- 1000
    for (ipow in 1:10) {
      set.seed(ipow)
      beta0 <- rep(uu[3], nrow(pattern))
      lambda0 <- sapply(1:numSignal, function(x) {
        exp(beta0 + rnorm(length(beta0), 0, itau/100))
      })
      
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:(numGene - numSignal), function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau/100))
      })
      newCt1 <- lapply(1:(numGene - numSignal), function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(rbind(do.call(rbind, newCt0), do.call(rbind, 
                                                                    newCt1)))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      write.csv(t(countdata), file = paste0("/Users/ke/Desktop/Practicum Project/made_by_myself/simulated_SE/sim_MOB_pattern", 
                                            ipt, "_fc", ifc, "_tau", itau, "_count_power", ipow, ".csv"), 
                row.names = T)
    }
  }
}

write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/simulation_ke/subset_simu/Rep11_MOB_info_spark.csv", row.names = T)


#### Model 1 Input simulated data into SPARK function:


## Analyze with SPARK
##-----------------
rm(list = ls())

library(SPARK)
itau = 35
ifc = 3
ipt = 2
ipow = 1

load(paste0("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_spark.rds")) 

info <- spark@location
tol_counts <- spark@lib_size
set.seed(123) # set seed for reproducibility
nor_count <- round(rnorm(length(tol_counts), mean(tol_counts), 10))


info <- read.csv("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/processed_data/Rep11_MOB_info_spark.csv", row.names = 1)

countdata <- t(read.csv(paste0("/Users/ke/Desktop/Practicum Project/made_by_myself/simulated_SE/sim_MOB_pattern", ipt, "_fc", ifc, 
                               "_tau", itau, "_count_power", ipow, ".csv"), row.names = 1))
spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], 
                           percentage = 0.1, min_total_counts = 10)

spark@lib_size <- nor_count

t1 <- proc.time()
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 1, verbose = F, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = F)
time_comp <- proc.time() - t1



#### Model 1 plot QQplot:

NumGene <- 1000

observed_pval <- spark@res_mtest$combined_pvalue
uni_dis <- seq(1, 1/NumGene, length.out = NumGene)

#take log
log_observed <- -log10(observed_pval) 
log_uni <- -log10(uni_dis)

qqplot(log_uni, log_observed, xlim=c(0,5), ylim=c(0,5), xlab="Expected -log10 p-value", ylab="Observed -log10 p-value", main="QQ plot of the observed -log10(P-value) from SPARK \n against the expected -log10(P-value)")
abline(0,1,col=2)



###############
##__Model 2__##
###############


#### Model 2 sum up cell counts as N:

##------------------------------------------------
## Data Generation
##------------------------------------------------
## Generate the countdata based on the info from the realdata and
## patterns summarized from the SpatialDE result

## Fold Change
##-----------------
rm(list = ls())

load(paste0("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_spark.rds")) 

info <- spark@location
tol_counts <- spark@lib_size
# random_number <- round(ceiling(rnorm(260, 29899.95, 12990.51)))

beta <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$coefficients
})

nb <- sapply(1:4, function(x) {
  log(x * exp(median(beta)))
})

tau2 <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$theta[2]
})


#### Model 2 Same as model 1 on generating total count:

set.seed(123) # set seed for reproducibility
nor_count <- round(rnorm(length(tol_counts), mean(tol_counts), 10))


#### Model 2 generate 11247 genes, including 1125 SE and 10122 non-SE genes, to make sure total number of genes is the same as the original real dataset. :

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35

for (ifc in 2:2) { #2:4
  newN <- nor_count 
  for (ipt in 2:2) { #2:4
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    
    # nrow(spark@counts) Totally 11274 genes
    numSignal <- 1125
    numGene <- 10122
    
    for (ipow in 1:1) {
      set.seed(ipow)
      beta0 <- rep(uu[3], nrow(pattern))
      lambda0 <- sapply(1:numSignal, function(x) {
        exp(beta0 + rnorm(length(beta0), 0, itau/100))
      })
      
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:(numGene - numSignal), function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau/100))
      })
      newCt1 <- lapply(1:(numGene - numSignal), function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(rbind(do.call(rbind, newCt0), do.call(rbind, 
                                                                    newCt1)))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      # Sum up for each spot (260)
      sum_count <- colSums(countdata)
      
    }
  }
}



#### Model 2 After sum up row counts, we simulated data:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35

for (ifc in 3:3) {
  newN <- sum_count
  for (ipt in 2:2) {
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    numSignal <- 100
    numGene <- 1000
    for (ipow in 1:10) {
      set.seed(ipow)
      beta0 <- rep(uu[3], nrow(pattern))
      lambda0 <- sapply(1:numSignal, function(x) {
        exp(beta0 + rnorm(length(beta0), 0, itau/100))
      })
      
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:(numGene - numSignal), function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau/100))
      })
      newCt1 <- lapply(1:(numGene - numSignal), function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(rbind(do.call(rbind, newCt0), do.call(rbind, 
                                                                    newCt1)))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      write.csv(t(countdata), file = paste0("/Users/ke/Desktop/Practicum Project/made_by_myself/simulated_NonSE/sim_MOB_pattern", 
                                            ipt, "_fc", ifc, "_tau", itau, "_count_power", ipow, ".csv"), 
                row.names = T)
    }
  }
}

write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/simulation_ke/subset_simu/Rep11_MOB_info_spark.csv", row.names = T)



#### Model 2 Input simulated data into SPARK function:

## Analyze with SPARK
##-----------------
rm(list=setdiff(ls(), "sum_count"))

library(SPARK)
itau = 35
ifc = 3
ipt = 2
ipow = 1

info <- read.csv("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/processed_data/Rep11_MOB_info_spark.csv", row.names = 1)

countdata <- t(read.csv(paste0("/Users/ke/Desktop/Practicum Project/made_by_myself/simulated_NonSE/sim_MOB_pattern", ipt, "_fc", ifc, "_tau", itau, "_count_power", ipow, ".csv"), row.names = 1))

spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], percentage = 0.1, min_total_counts = 10)

spark@lib_size <- sum_count

t1 <- proc.time()
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 1, verbose = F, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = F)
time_comp <- proc.time() - t1



#### Model 2 plot QQplot:

NumGene <- 1000

observed_pval <- spark@res_mtest$combined_pvalue
uni_dis <- seq(1, 1/NumGene, length.out = NumGene)

#take log
log_observed <- -log10(observed_pval) 
log_uni <- -log10(uni_dis)

qqplot(log_uni, log_observed, xlim=c(0,5), ylim=c(0,5), xlab="Expected -log10 p-value", ylab="Observed -log10 p-value", main="QQ plot of the observed -log10(P-value) from SPARK \n against the expected -log10(P-value)")
abline(0,1,col=2)




