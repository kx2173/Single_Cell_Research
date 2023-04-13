
##########################
# Power Plot of Mob Data #
##########################



library(SPARK)
library(ggplot2)


##------------------------------------------------
## Visualization of the pattern summarized from SpatialDE result
##------------------------------------------------

rm(list = ls())
source("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/funcs/funcs.R")

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

tmplist <- datlist
for (i in 1:5) {
  tmplist[[i]][, 2] <- relative_func(datlist[[i]][, 2])
}

patterns = c("I", "II", "III")
## three major pattern were used for simulation
df <- setNames(cbind.data.frame(tmplist[[1]][, 1], do.call(cbind, sapply(tmplist[c(5, 
                                                                                   1, 4)], "[", 2))), c("xy", paste0("Pattern ", patterns)))
pp <- lapply(1:3, function(x) {
  pattern_plot2(df, x, xy = F, main = T, titlesize = 1.5)
})

#pattern 2,3,4
grid.arrange(grobs = pp, ncol = 3)



#### Extract coefficients from SPARK object:

##------------------------------------------------
## Data Generation
##------------------------------------------------
## Generate the countdata based on the info from the realdata and
## patterns summarized from the SpatialDE result

## Fold Change
##-----------------
rm(list = ls())

load(paste0("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_spark.rds")) 

# Generate info with 260 obs and 3 variables (x, y, total_counts)
info <- spark@location # 260 spots with x-axis and y-axis
info$total_counts <- spark@lib_size # Counts for 260 spots

beta <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$coefficients #333 gene coefficient
})

nb <- sapply(1:4, function(x) {
  log(x * exp(median(beta)))
})

tau2 <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$theta[2]
})



#### Calculate SE proportion:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35
pvalue_SE <- data.frame()
alpha <- seq(0, .1, .0005)
SE_comb <- data.frame()

for (ifc in 2:4) {
  newN <- info$total_counts
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
      
      spark@lib_size <- info$total_counts
      
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, num_core = 1, verbose = F, fit.maxiter = 500) #returns tau for 10,000 gene
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue
      
      SE_comb <- rbind(SE_comb, comb_pvalue)
      
      rm(countdata, spark, spark_p, comb_pvalue)
    }
  }
}

write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)



#### Calculate SE power:

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



#### Calculate Non-SE proportion:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35
pvalue_NonSE <- data.frame()
alpha <- seq(0, .1, .0005)
NonSE_comb <- data.frame()

for (ifc in 2:4) {
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
      
      spark@lib_size <- info$total_counts
      
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, num_core = 1, verbose = F, fit.maxiter = 500) #returns tau for 10,000 gene
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue
      
      NonSE_comb <- rbind(NonSE_comb, comb_pvalue)
      
      rm(countdata, spark, spark_p, comb_pvalue)
    }
  }
}

write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)



#### Calculat FDR:

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



#### Draw Power plot:

plot(test_fdr, test_power, type = "l", col="pink", lwd=3,
     xlim=c(0,0.1), ylim=c(0,1),
     xlab="FDR", ylab="Power", main="Power plot at a range of FDR")










