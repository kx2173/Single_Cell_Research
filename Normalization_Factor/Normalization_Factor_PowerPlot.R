#######################################
# Power Plot for Normalization Factor #
#######################################


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


#### Model 1 Calculate SE proportion:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35
pvalue_SE <- data.frame()
alpha <- seq(0, .1, .0005)
SE_comb <- data.frame()

for (ifc in 2:2) {
  newN <- nor_count
  for (ipt in 3:3) {
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    numSignal <- 100
    numGene <- 1000
    for (ipow in 1:5) {
      set.seed(ipow)
      
      beta0 <- uu[grp]
      lambda0 <- sapply(1:numSignal, function(x) {
        exp(beta0 + rnorm(length(beta0), 0, itau/100))
      })
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      countdata <- data.frame(do.call(rbind, newCt0))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], percentage = 0.1, min_total_counts = 10)
      
      spark@lib_size <- nor_count
      
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, num_core = 1, verbose = F, fit.maxiter = 500) #returns tau for 10,000 gene
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue
      
      SE_comb <- rbind(SE_comb, comb_pvalue)
      
      rm(countdata, spark, spark_p, comb_pvalue)
    }
  }
}

# write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)



#### Model 1 Calculate SE power:

test_se <- data.frame()

SE_comb_new <- t(SE_comb)

for (j in 1:ncol(SE_comb_new)) {
  for (i in 1:length(alpha)) {
    test_se[i,j] <- sum(SE_comb_new[,j] < alpha[i]) / numSignal
    i <- i+1
  }
  j <- j+1
}

test_power <- rowMeans(test_se)


#### Model 1 Calculate Non-SE proportion:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35
pvalue_NonSE <- data.frame()
alpha <- seq(0, .1, .0005)
NonSE_comb <- data.frame()

for (ifc in 2:2) {
  for (ipt in 3:3) {
    for (ipow in 1:5) {
      set.seed(ipow)
      
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:(numGene - numSignal), function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau/100))
      })
      newCt1 <- lapply(1:(numGene - numSignal), function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(do.call(rbind, newCt1))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], percentage = 0.1, min_total_counts = 10)
      
      spark@lib_size <- nor_count
      
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, num_core = 1, verbose = F, fit.maxiter = 500) #returns tau for 10,000 gene
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue
      
      NonSE_comb <- rbind(NonSE_comb, comb_pvalue)
      
      rm(countdata, spark, spark_p, comb_pvalue)
    }
  }
}

# write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)



#### Model 1 Non SE calculation:

test_p <- data.frame()

NonSE_comb_new <- t(NonSE_comb)

for (j in 1:ncol(NonSE_comb_new)) {
  for (i in 1:length(alpha)) {
    test_p[i,j] <- sum(NonSE_comb_new[,j] < alpha[i]) / (numGene - numSignal)
    i <- i+1
  }
  j <- j+1
}

test_fdr <- rowMeans(test_p)


#### Model 1 Plot power plot:

plot(test_fdr, test_power, type = "l", col="pink", lwd=3,
     xlim=c(0,0.1), ylim=c(0,1),
     xlab="FDR", ylab="Power", main="Power plot at a range of FDR")






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
rm(list=setdiff(ls(), c("test_fdr", "test_power")))

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



#### Model 2 generate 11247 genes, including 1125 SE and 10122 non-SE genes, to make sure total number of genes is the same as the original real dataset :

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35

for (ifc in 2:2) { #2:4
  newN <- nor_count 
  for (ipt in 3:3) { #2:4
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



#### Model 2 Calculate SE proportion:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35
pvalue_SE <- data.frame()
alpha <- seq(0, .1, .0005)
SE_comb <- data.frame()

for (ifc in 2:2) {
  newN <- sum_count
  for (ipt in 3:3) {
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    numSignal <- 100
    numGene <- 1000
    for (ipow in 1:5) {
      set.seed(ipow)
      
      beta0 <- uu[grp]
      lambda0 <- sapply(1:numSignal, function(x) {
        exp(beta0 + rnorm(length(beta0), 0, itau/100))
      })
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      countdata <- data.frame(do.call(rbind, newCt0))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], percentage = 0.1, min_total_counts = 10)
      
      spark@lib_size <-sum_count
      
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, num_core = 1, verbose = F, fit.maxiter = 500) #returns tau for 10,000 gene
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue

      SE_comb <- rbind(SE_comb, comb_pvalue)
      
      rm(countdata, spark, spark_p, comb_pvalue)
    }
  }
}

write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)



#### Model 2 Calculate SE power:

test_se <- data.frame()

SE_comb_new <- t(SE_comb)

for (j in 1:ncol(SE_comb_new)) {
  for (i in 1:length(alpha)) {
    test_se[i,j] <- sum(SE_comb_new[,j] < alpha[i]) / numSignal
    i <- i+1
  }
  j <- j+1
}

test_power2 <- rowMeans(test_se)



#### Model 2 Calculate Non-SE proportion:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35
pvalue_NonSE <- data.frame()
alpha <- seq(0, .1, .0005)
NonSE_comb <- data.frame()

for (ifc in 2:2) {
  for (ipt in 3:3) {
    for (ipow in 1:5) {
      set.seed(ipow)
      
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:(numGene - numSignal), function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau/100))
      })
      newCt1 <- lapply(1:(numGene - numSignal), function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(do.call(rbind, newCt1))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], percentage = 0.1, min_total_counts = 10)
      
      spark@lib_size <- sum_count
      
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, num_core = 1, verbose = F, fit.maxiter = 500) #returns tau for 10,000 gene
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue
      
      NonSE_comb <- rbind(NonSE_comb, comb_pvalue)
      
      rm(countdata, spark, spark_p, comb_pvalue)
    }
  }
}

write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)




#### Model 2 Non SE calculation:

test_p <- data.frame()

NonSE_comb_new <- t(NonSE_comb)

for (j in 1:ncol(NonSE_comb_new)) {
  for (i in 1:length(alpha)) {
    test_p[i,j] <- sum(NonSE_comb_new[,j] < alpha[i]) / (numGene - numSignal)
    i <- i+1
  }
  j <- j+1
}

test_fdr2 <- rowMeans(test_p)



#### Model 2 Plot Power plot:

plot(test_fdr2, test_power2, type = "l", col="pink", lwd=3,
     xlim=c(0,0.1), ylim=c(0,1),
     xlab="FDR", ylab="Power", main="Power plot at a range of FDR")







