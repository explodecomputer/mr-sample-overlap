# Load required libraries
# Need devtools to be installed
# install.packages("devtools")
# library(devtools)
# install_github("MRCIEU/TwoSampleMR")

library(simulateGP)
library(TwoSampleMR)
#library(progress)


# Function for running the algorithm on bc4
main <- function()
{
  param <- expand.grid(
    nid = 20000,
    nsnp = 20,
    maf = 0.3,
    nsim = 1:10000,
    xy = c(0, 0.2),
    ux = 1,
    uy = c(0, 0.6, 1, 2),
    gx = seq(0.04, 0.24, by=0.04),
    offset = seq(0, 10000, length.out=11)
  )
  
  args <- commandArgs(T)
  splitsize <- as.numeric(args[1])
  chunk <- as.numeric(args[2])
  set.seed(chunk)
  
  first <- (chunk-1) * splitsize + 1
  last <- min(chunk * splitsize, nrow(param))
  param <- param[first:last,]
  
  message("Running ", nrow(param), " simulations")
  message(first, " to ", last)
  
  out <- list()
  #	pb <- progress_bar$new(total=nrow(param))
  for(i in 1:nrow(param))
  {
    #		pb$tick()
    message(i)
    out[[i]] <- suppressMessages(runsim(param[i,]))
  }
  out <- bind_rows(out)
  dir.create("../results/sim5", recursive=TRUE, showWarnings=FALSE)
  save(out, file=paste0("../results/sim5/out", chunk, ".rdata"))
}

# Function for running the simulations
runsim <- function(param)
{
  p2 <- param
  g <- make_geno(param$nid, param$nsnp, param$maf)
  u <- rnorm(param$nid)
  x <- cbind(g,u) %*% c(rep(param$gx, param$nsnp), param$ux) + rnorm(param$nid)
  y <- param$xy * x + u * param$uy + rnorm(param$nid)
  
  # Effect of G on X in dataset D
  exp_d <- gwas(x[1:(param$nid/2)], g[1:(param$nid/2),])
  # Effect of G on X in replication R
  exp_r <- gwas(x[(param$nid/2+1):param$nid], g[(param$nid/2+1):param$nid,])
  # Effect of G on Y
  out <- gwas(y[1:(param$nid/2)+param$offset], g[1:(param$nid/2)+param$offset,])
  
  # perform the MR analysis using D
  dat <- data.frame(
    SNP = 1:nrow(exp_d),
    beta.exposure = exp_d$bhat,
    beta.outcome = out$bhat,
    se.exposure = exp_d$se,
    se.outcome = out$se,
    pval.exposure = exp_d$pval,
    pval.outcome = out$pval,
    fval.exposure = exp_d$fval,
    id.exposure=1,
    id.outcome=2,
    exposure=1,
    outcome=1,
    mr_keep=TRUE
  )
  res <- suppressMessages(mr(dat, method=c("mr_wald_ratio", "mr_ivw")))
  
  # Why?
  param$b <- res$b
  param$se <- res$se
  param$nsnp_inc <- nrow(dat)
  param$pval <- res$pval
  # This causing Error: object of type 'builtin' is not subsettable
  param$mean_f <- mean(exp_d$fval, na.rm=TRUE)
  
  param$obs <- lm(y ~ x)$coef[2]
  param$what <- "all"
  
  # filter by most sifnificant SNPs (in mr or dat?)
  dat_s <- subset(dat, pval.exposure < 5e-8)
  
  # SNPs from the exposure replication dataset that gave significant results
  exp_r_s <- subset(exp_r, snp %in% dat_s$SNP)
  
  # SNPs from the outcome dataset that correspond to the significant results
  out_s <- subset(out, snp %in% dat_s$SNP)
  
  # re-do the MR analysis for the significant SNPs, if there are any
  if (nrow(dat_s) > 0)
  {
    dat <- data.frame(
      SNP = 1:nrow(exp_r_s),
      beta.exposure = exp_r_s$bhat,
      beta.outcome = out_s$bhat,
      se.exposure = exp_r_s$se,
      se.outcome = out_s$se,
      pval.exposure = exp_r_s$pval,
      pval.outcome = out_s$pval,
      fval.exposure = exp_r_s$fval,
      id.exposure=1,
      id.outcome=2,
      exposure=1,
      outcome=1,
      mr_keep=TRUE
    )
    res <- suppressMessages(mr(dat, method=c("mr_wald_ratio", "mr_ivw")))
  }
  
  
  
  if(nrow(dat) > 0)
  {	
    res2 <- suppressMessages(mr(dat, method="mr_wald_ratio", "mr_ivw"))
    p2$b <- res$b
    p2$se <- res$se
    p2$nsnp_inc <- nrow(dat)
    p2$pval <- res$pval
    p2$mean_f <- mean(dat$fval.exposure, na.rm=TRUE)
    p2$obs <- lm(y ~ x)$coef[2]
  }
  p2$what <- "sig"
  out <- bind_rows(param, p2)
  return(out)
}


param <- expand.grid(
  nid = 20000,
  nsnp = 20,
  maf = 0.3,
  nsim = 1:10000,
  xy = c(0, 0.2),
  ux = 1,
  uy = c(0, 0.6, 1, 2),
  gx = seq(0.04, 0.24, by=0.04),
  offset = seq(0, 10000, length.out=11)
)

runsim(param[100000,])

dim(param)
main()

