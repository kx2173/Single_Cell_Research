library(SPARK)
library(foreach)
library(doParallel)



# Visualization of pattern from SpatialDE result:

##------------------------------------------------
## Visualization of the pattern summarized from SpatialDE result
##------------------------------------------------

# rm(list = ls())
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

#pattern 2,3,4 -> ipt=1,2,3
# grid.arrange(grobs = pp, ncol = 3)

pp[2]




# Model 1 generate 100 SE and 900 non-SE genes:
rm(list = ls())
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
        exp(beta0 + rnorm(length(beta0), 0, itau))
      })
      
      newCt0 <- lapply(1:numSignal, function(x) {
        rpois(length(lambda0[, x]), newN * lambda0[, x])
      })
      
      beta1 <- rep(uu[3], nrow(pattern))
      lambda1 <- sapply(1:(numGene - numSignal), function(x) {
        exp(beta1 + rnorm(length(beta1), 0, itau))
      })
      newCt1 <- lapply(1:(numGene - numSignal), function(x) {
        rpois(length(lambda1[, x]), newN * lambda1[, x])
      })
      
      countdata <- data.frame(rbind(do.call(rbind, newCt0), do.call(rbind, 
                                                                    newCt1)))
      rownames(countdata) <- paste0("gene", 1:nrow(countdata))
      colnames(countdata) <- pattern[, 1]
    }
  }
}

# workdir <- "/Users/sw3763/Documents/Research/scRNA-seq/MS project/Ke Xu/04-11-2023/"
# write.csv(info, file = paste0(workdir, "Rep11_MOB_info_spark.csv"), row.names = T)




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



spark0 <- CreateSPARKObject(counts = countdata, location = info[, 1:2], 
                            percentage = 0.1, min_total_counts = 10)

spark0@lib_size <- nor_count

# True normalization factor
t1 <- proc.time()
spark <- spark.vc(spark0, covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 1, verbose = F, fit.maxiter = 500)
spark <- spark.test(spark, check_positive = T, verbose = F)
time_comp <- proc.time() - t1

# Simple average normalization factor
t1 <- proc.time()
spark3 <- spark.vc(spark0, covariates = NULL, lib_size = NULL, 
                   num_core = 1, verbose = F, fit.maxiter = 500)
spark3 <- spark.test(spark3, check_positive = T, verbose = F)
time_comp3 <- proc.time() - t1

# Moment based normalization factor

t1 <- proc.time()
spark2 <- spark.vc2(spark0, covariates = NULL, lib_size = spark@lib_size, 
                    num_core = 1, verbose = F, fit.maxiter = 500)
spark2 <- spark.test(spark2, check_positive = T, verbose = F)
time_comp2 <- proc.time() - t1



## QQplot:
NumGene <- 1000

observed_pval <- spark@res_mtest$combined_pvalue
observed_pval2 <- spark2@res_mtest$combined_pvalue
observed_pval3 <- spark3@res_mtest$combined_pvalue
uni_dis <- seq(1, 1/NumGene, length.out = NumGene)

#take log
log_observed <- -log10(observed_pval) 
log_observed2 <- -log10(observed_pval2)
log_observed3 <- -log10(observed_pval3)
log_uni <- -log10(uni_dis)

qqplot(log_uni, log_observed, xlim=c(0,3), ylim=c(0,3), pch=17, col=2, 
       xlab="Expected -log10 p-value", ylab="Observed -log10 p-value", 
       main="QQ plot of the observed -log10(P-value) from SPARK \n against the expected -log10(P-value)", frame.plot = FALSE)
axis(side = 1)
par(new=TRUE)
qqplot(log_uni, log_observed2, xlim=c(0,3), ylim=c(0,3), pch=18, col=4, 
       xlab="", ylab="", frame.plot = FALSE)
par(new=TRUE)
qqplot(log_uni, log_observed3, xlim=c(0,3), ylim=c(0,3), pch=20, col='purple', 
       xlab="", ylab="", frame.plot = FALSE)
abline(0,1,col=8)
legend('topleft',c('True','Estimate','Simple avg'),pch=c(17,18,20), 
       col = c(2,4,'purple'), bty = "n")



Interpretation: 
  Quantile–quantile plot of the observed −log10(P) from different methods against the expected −log10(P) under the null condition for the first set of null simulations. Simulations were performed under moderate noise ($\tau_2 = 0.35$). Compared normalization factor include 

- True (red triangle)

- Estimated (blue diamond)

- Simple average (purple dot)

Findings: 
  Under the null condition, Estimated normalization factor produced well-calibrated P values at the transcriptome-wide significance levels. Simple average normalization factor produced test statistics yielded slightly conservative P values. 

Explanation:
  The failure of using simple average normalization factor in controlling type I errors is presumably due to its ignorance of the variation from residual noise.

<!-- The P value calibration results under the null condition for different methods were consistent across simulation settings and across a range of noise variance levels. -->
  
  <!-- A representative expression pattern for a null gene is embedded inside the panel. -->
  
  ## Power plot (need to compare 3 models)
