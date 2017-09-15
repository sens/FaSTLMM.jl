##################################################################
## R/lmmlite results and comparison with Julia/FaSTLMM and Pylmm
##################################################################
## Builds on code by Karl Broman
##################################################################

# load library and attach data
library(lmmlite)
data(recla)
n_phe <- ncol(recla$pheno)
lmmlite <- data.frame(index=rep(NA,2*n_phe),
                      method=rep("", 2*n_phe),
                      hsq=rep(NA,2*n_phe),
                      intercept=rep(NA,2*n_phe),
                      sex=rep(NA,2*n_phe),
                      sigmasq=rep(NA,2*n_phe),
                      loglik=rep(NA,2*n_phe), stringsAsFactors=FALSE)

## go through all the phenotypes and get results

all_e <- vector("list", n_phe)
for(i in 1:n_phe) {
    keep <- !is.na(recla$pheno[,i])
    y <- recla$pheno[keep,i,drop=FALSE]
    x <- recla$covar[keep,]
    k <- recla$kinship[keep,keep]

    all_e[[i]] <- e <- eigen_rotation(k, y, x)
    res <- fitLMM(e$Kva, e$y, e$X, reml=TRUE)
    lmmlite[i*2-1,1] <- i
    lmmlite[i*2-1,2] <- "reml"
    lmmlite[i*2-1,-(1:2)] <- c(res$hsq, res$beta, res$sigmasq, res$loglik)

    res <- fitLMM(e$Kva, e$y, e$X, reml=FALSE)
    lmmlite[i*2,1] <- i
    lmmlite[i*2,2] <- "ml"
    lmmlite[i*2,-(1:2)] <- c(res$hsq, res$beta, res$sigmasq, res$loglik)
}

## write results to file
write.csv(lmmlite,"lmmlite_results.csv")

#################################################################
## performance comparison
#################################################################

library(microbenchmark)
library(broman)

time_full2 <- microbenchmark(R={
    for(i in 1:ncol(recla$pheno)) {
        keep <- !is.na(recla$pheno[,i])
        eigenrot <- eigen_rotation(recla$kinship[keep,keep],
                                   recla$pheno[keep,i],
                                   recla$covar[keep,], use_cpp=FALSE)
        res <- fitLMM(eigenrot$Kva, eigenrot$y, eigenrot$X, reml=TRUE,
                      use_cpp=FALSE)
    }},
    cpp={
        for(i in 1:ncol(recla$pheno)) {
            keep <- !is.na(recla$pheno[,i])
            eigenrot <- eigen_rotation(recla$kinship[keep,keep],
                                       recla$pheno[keep,i],
                                       recla$covar[keep,], use_cpp=TRUE)
            res <- fitLMM(eigenrot$Kva, eigenrot$y, eigenrot$X, reml=TRUE,
                          use_cpp=TRUE)
        }}, times=100)
print(time_full2, digits=4)

