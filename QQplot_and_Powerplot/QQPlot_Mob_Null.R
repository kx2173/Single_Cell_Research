##################
# QQ Plot of Mob #
##################


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


grid.arrange(grobs = pp, ncol = 3)




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
info$total_counts <- spark@lib_size

beta <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$coefficients
})

nb <- sapply(1:4, function(x) {
  log(x * exp(median(beta)))
})

tau2 <- sapply(1:length(spark@res_vc), function(x) {
  spark@res_vc[[x]]$theta[2]
})




#### Simulated data under NULL conditions:

load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

itau = 35

for (ifc in 2:4) {
  newN <- info$total_counts
  for (ipt in 2:4) {
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    numGene <- 1000
    for (ipow in 1:3) {
      set.seed(ipow)
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:numGene, function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau/100))
      })
      newCt1 <- lapply(1:numGene, function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(do.call(rbind, newCt1))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      write.csv(t(countdata), file = paste0("/Users/ke/Desktop/Practicum Project/made_by_myself/simulated_SE/sim_MOB_pattern", 
                                            ipt, "_fc", ifc, "_tau", itau, "_count_power", ipow, ".csv"), 
                row.names = T)
    }
  }
}

write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/simulation_ke/subset_simu/Rep11_MOB_info_spark.csv", row.names = T)



#### Input simulated data into SPARK function:

## Analyze with SPARK
##-----------------
rm(list = ls())

library(SPARK)
itau = 35
ifc = 3
ipt = 2
ipow = 1

info <- read.csv("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/processed_data/Rep11_MOB_info_spark.csv", row.names = 1)

countdata <- t(read.csv(paste0("/Users/ke/Desktop/Practicum Project/made_by_myself/simulated_SE/sim_MOB_pattern", ipt, "_fc", ifc, 
                               "_tau", itau, "_count_power", ipow, ".csv"), row.names = 1))
spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], 
                           percentage = 0.1, min_total_counts = 10)

spark@lib_size <- info$total_counts
t1 <- proc.time()
spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 1, verbose = F, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = F)
time_comp <- proc.time() - t1




#### Draw QQplot:

NumGene <- 1000

observed_pval <- spark@res_mtest$combined_pvalue
uni_dis <- seq(1, 1/NumGene, length.out = NumGene)

#take log
log_observed <- -log10(observed_pval) 
log_uni <- -log10(uni_dis)

qqplot(log_uni, log_observed, xlim=c(0,5), ylim=c(0,5), xlab="Expected -log10 p-value", ylab="Observed -log10 p-value", main="QQ plot of the observed -log10(P-value) from SPARK \n against the expected -log10(P-value)")
abline(0,1,col=2)






