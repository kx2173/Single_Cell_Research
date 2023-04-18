library(SPARK)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)


# spark.vc2:
spark.vc2 <- function(object,
                      covariates=NULL,
                      lib_size=NULL,
                      fit.maxiter=500,
                      fit.tol=1e-5,
                      fit.model="poisson",
                      num_core=1,
                      verbose=FALSE) {
  ## extract the data from the slot of object, createobject() function goes first
  if(length(object@counts) == 0) {
    stop("object@counts has not been set. Run CreateSPARKObject() first and then retry.")
  }# end fi
  
  ## load necessary R packages
  #require("CompQuadForm")
  #require("doParallel")
  #require("foreach")
  #require("Matrix")
  
  ## set number of core in run
  if(num_core > 1){
    if(num_core>detectCores()){warning("SPARK.VC:: the number of cores you're setting is larger than detected cores!");
      num_core = detectCores()}
  }# end fi
  
  ## store the multiple cores
  object@num_core <- num_core
  ## covariates, i.e., confounding or batch effects
  if(is.null(covariates)){
    num_cov <- 0
  }else{
    # remove the intercept if added by user, later intercept will add automatically
    if(length(unique(covariates[,1])) == 1){
      covariates <- covariates[, -1]
    }# end fi
    covariates <- as.matrix(covariates)
    num_cov <- ncol(covariates)
  }# end fi
  
  
  ## number of genes and cells
  num_gene <- nrow(object@counts)
  num_cell <- ncol(object@counts)
  
  cat(paste("## ===== SPARK INPUT INFORMATION ==== \n"))
  cat(paste("## number of total samples: ", num_cell,"\n"))
  cat(paste("## number of total features: ", num_gene,"\n"))
  cat(paste("## number of adjusted covariates: ", num_cov,"\n"))
  
  object@fit_model <- fit.model
  ## main functions
  if(fit.model=="poisson"){
    #*************************************************#
    # Count-Based Spatial Model Under The Null #
    #*************************************************#
    cat("# fitting count-based spatial model under the null hypothesis ... \n")
    if(is.null(lib_size)){
      lib_size <- apply(object@counts, 2, sum)
    }else{# already exists
      lib_size <- as.numeric(lib_size)
    }# end fi
    object@lib_size <- lib_size
    ig <- 1
    
    
    lib_size.new <- lib_size
    counts.new <- object@counts
    for(Iter in 1:2) {
      lib_size.old <- lib_size.new
      counts.old <- counts.new
      counts.new <- foreach(ig = 1:num_gene,.combine = 'rbind')%do%{
        model0 <- try(glm(formula = as.numeric(counts.old[ig,]) ~ 1 +
                            offset(log(lib_size.old)), family = quasipoisson(link="log")))
        model0$fitted.values
      }
      lib_size.new <- apply(counts.new, 2, sum)
      cat('\nlib_size_diff=',max(abs(lib_size.new-lib_size.old)))
    }
    lib_size <- lib_size.new
    
    res_vc <- foreach(ig = 1:num_gene)%do%{
      if(num_cov==0){
        model0 <- try(glm(formula = as.numeric(object@counts[ig,]) ~ 1 +
                            offset(log(lib_size)), family = poisson(link="log")))
        idx <- match(rownames(model.frame(formula = as.numeric(object@counts[ig,]) ~ 1 +
                                            offset(log(lib_size)), na.action = na.omit)),
                     rownames(model.frame(formula = as.numeric(object@counts[ig,]) ~ 1 +
                                            offset(log(lib_size)), na.action = na.pass)))
      }else{
        model0 <- try(glm(formula = as.numeric(object@counts[ig,]) ~ covariates +
                            offset(log(lib_size)), family = poisson(link="log")))
        idx <- match(rownames(model.frame(formula = object@counts[ig,] ~ covariates +
                                            offset(log(lib_size)), na.action = na.omit)),
                     rownames(model.frame(formula = object@counts[ig,] ~ covariates +
                                            offset(log(lib_size)), na.action = na.pass)))
      }# end fi
      
      ## model to fit
      if(verbose) {cat(paste("NO. Gene = ",ig,"\n"))}
      t1 <- system.time(model1 <- try(spark.fit(model0, maxiter = fit.maxiter,
                                                tol = fit.tol, verbose=verbose)))
      
      ## store results
      if((class(model1) != "try-error")&&(!any(is.na(model1$Y)))){
        if(verbose){cat(paste("SPARK.CV::tau = ", model1$theta,"\n"))}
        model1$analy_indx <- idx # cell index to run for each gene
        model1$times <- t1[3] # running time for each gene
      }else{
        model1$converged <- FALSE
      }# end fi
      #######
      model1
      #res_vc[[ig]] <- model1
      ########
    }# end for ig, parallel
    #close(pb)
    #stopCluster(cl)
    
    names(res_vc) <- rownames(object@counts)
    object@res_vc <- res_vc
  }else if(fit.model=="gaussian"){
    #************************************************************#
    # Normalized Count-Based Spatial Model Under The Null #
    #************************************************************#
    cat("# fitting normalized count-based spatial model under the null hypothesis ... \n")
    ## check normalizd counts, vst normalization method
    if(is.null(object@scaled_counts)){
      object@scaled_counts <- NormalizeVST(object@counts)
    }# end fi
    
    ## covariates are required in the downstream steps
    ## if null add a column vector with 1
    if(is.null(covariates)){
      covariates <- as.matrix(rep(1, num_cell))
    }# end fi
    
    object@res_vc <- list(covariates=as.matrix(covariates))
  }# end fi
  
  # return results
  return(object)
}# end function




# Model 1 Calculate SE proportion:
# rm(list = ls())
datapath <- '/Users/ke/Desktop/'
load(paste0(datapath,"Rep11_MOB_spark.rds")) 

# nor_count <- round(rnorm(length(tol_counts), mean(tol_counts), 10))
info <- spark@location
tol_counts <- spark@lib_size
load(paste0(datapath,"/MOB_Pattern_SpatialDE.rds"),verbose=T)

set.seed(123) # set seed for reproducibility
itau <- 0.35
ifc = 3
ipt = 2
ipow = 1
nb <- c(-10.2,-9.5,-9.1,-8.8)
nor_count <- runif(260,10000,30000)

pvalue_SE <- data.frame()
alpha <- seq(0, .1, .0005)
SE_comb <- SE_comb3 <- SE_comb2 <- data.frame()

for (ifc in 2:2) {
  newN <- nor_count
  for (ipt in 3:3) {
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    uu <- c(nb[1], nb[ifc], nb[1])
    numSignal <- 100
    numGene <- 1000
    for (ipow in 1:1) {
      set.seed(ipow)
      
      beta0 <- uu[grp]
      lambda0 <- sapply(1:numSignal, function(x) {
        exp(beta0 + rnorm(length(beta0), 0, itau))
      })
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      countdata <- data.frame(do.call(rbind, newCt0))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
      
      spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], 
                                 percentage = 0.1, min_total_counts = 10)
      
      spark@lib_size <- nor_count
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                        num_core = 1, verbose = F, fit.maxiter = 500) 
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue
      SE_comb <- rbind(SE_comb, comb_pvalue)
      
      spark3 <- spark.vc(spark, covariates = NULL, lib_size = NULL, 
                         num_core = 1, verbose = F, fit.maxiter = 500)
      spark_p3 <- spark.test(spark3, check_positive = T, verbose = F)
      comb_pvalue3 <- spark_p3@res_mtest$combined_pvalue
      SE_comb3 <- rbind(SE_comb3, comb_pvalue3)
      
      spark2 <- spark.vc2(spark, covariates = NULL, lib_size = spark@lib_size, 
                          num_core = 1, verbose = F, fit.maxiter = 500)
      spark_p2 <- spark.test(spark2, check_positive = T, verbose = F)
      comb_pvalue2 <- spark_p2@res_mtest$combined_pvalue
      SE_comb2 <- rbind(SE_comb2, comb_pvalue2)
      
      rm(countdata, spark, spark_p, comb_pvalue, 
         spark3, spark_p3, comb_pvalue3, spark2, spark_p2, comb_pvalue2)
      
    }
  }
}

# write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)



# Model 1 Calculate SE power:
test_se <- test_se2 <- test_se3 <- data.frame()

SE_comb_new <- t(SE_comb)
SE_comb_new2 <- t(SE_comb2)
SE_comb_new3 <- t(SE_comb3)

for (j in 1:ncol(SE_comb_new)) {
  for (i in 1:length(alpha)) {
    test_se[i,j] <- sum(SE_comb_new[,j] < alpha[i]) / numSignal
    test_se2[i,j] <- sum(SE_comb_new2[,j] < alpha[i]) / numSignal
    test_se3[i,j] <- sum(SE_comb_new3[,j] < alpha[i]) / numSignal
    i <- i+1
  }
  j <- j+1
}

test_power <- rowMeans(test_se)
test_power2 <- rowMeans(test_se2)
test_power3 <- rowMeans(test_se3)



# Model 1 Calculate Non-SE proportion:
load("/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/output/MOB_Pattern_SpatialDE.rds")

datapath <- '/Users/ke/Desktop/R_SPARK/utilities/'
load(paste0(datapath,"Rep11_MOB_spark.rds")) 

itau = 35
pvalue_NonSE <- data.frame()
alpha <- seq(0, .1, .0005)
NonSE_comb <- NonSE_comb2 <- NonSE_comb3 <- data.frame()

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
      
      spark <- CreateSPARKObject(counts = countdata, location = info[, 1:2], 
                                 percentage = 0.1, min_total_counts = 10)
      
      spark@lib_size <- nor_count
      spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                        num_core = 1, verbose = F, fit.maxiter = 500) 
      
      spark_p <- spark.test(spark, check_positive = T, verbose = F)
      comb_pvalue <- spark_p@res_mtest$combined_pvalue
      NonSE_comb <- rbind(NonSE_comb, comb_pvalue)
      
      spark3 <- spark.vc(spark, covariates = NULL, lib_size = NULL, 
                         num_core = 1, verbose = F, fit.maxiter = 500)
      spark_p3 <- spark.test(spark3, check_positive = T, verbose = F)
      comb_pvalue3 <- spark_p3@res_mtest$combined_pvalue
      NonSE_comb3 <- rbind(NonSE_comb3, comb_pvalue3)
      
      spark2 <- spark.vc2(spark, covariates = NULL, lib_size = spark@lib_size, 
                          num_core = 1, verbose = F, fit.maxiter = 500)
      spark_p2 <- spark.test(spark2, check_positive = T, verbose = F)
      comb_pvalue2 <- spark_p2@res_mtest$combined_pvalue
      NonSE_comb2 <- rbind(NonSE_comb2, comb_pvalue2)
      
      rm(countdata, spark, spark_p, comb_pvalue, 
         spark3, spark_p3, comb_pvalue3, spark2, spark_p2, comb_pvalue2)
      
    }
  }
}

# write.csv(info, file = "/Users/ke/Desktop/Practicum Project/project_spark/SPARK-Analysis-master/test_ke/Rep11_MOB_info_spark.csv", row.names = T)



# Model 1 Calculate Non-SE power:
test_non <- test_non2 <- test_non3 <- data.frame()

SE_Non_new <- t(NonSE_comb)
SE_Non_new2 <- t(NonSE_comb2)
SE_Non_new3 <- t(NonSE_comb3)

for (j in 1:ncol(SE_comb_new)) {
  for (i in 1:length(alpha)) {
    test_non[i,j] <- sum(SE_Non_new[,j] < alpha[i]) / (numGene - numSignal)
    test_non2[i,j] <- sum(SE_Non_new2[,j] < alpha[i]) / (numGene - numSignal)
    test_non3[i,j] <- sum(SE_Non_new3[,j] < alpha[i]) / (numGene - numSignal)
    i <- i+1
  }
  j <- j+1
}

test_fdr <- rowMeans(test_non)
test_fdr2 <- rowMeans(test_non2)
test_fdr3 <- rowMeans(test_non3)



# Combined plot:
qqplot(test_fdr, test_power, xlim=c(0,0.1), ylim=c(0,1), pch=17, col=2, 
       xlab="FDR", ylab="Power", 
       main="Power plot of three models with different factors \n under the null condition", frame.plot = FALSE)
axis(side = 1)
par(new=TRUE)
qqplot(test_fdr2, test_power2, xlim=c(0,0.1), ylim=c(0,1), pch=18, col=4, 
       xlab="", ylab="", frame.plot = FALSE)
par(new=TRUE)
qqplot(test_fdr3, test_power3, xlim=c(0,0.1), ylim=c(0,1), pch=20, col='purple', 
       xlab="", ylab="", frame.plot = FALSE)
legend('topleft',c('Model 1','Model 2','Model 3'),pch=c(17,18,20), 
       col = c(2,4,'purple'), bty = "n")


