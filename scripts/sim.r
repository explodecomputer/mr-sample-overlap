# Load required libraries
# Need devtools to be installed
# install.packages("devtools")
# library(devtools)
# install_github("MRCIEU/TwoSampleMR")

library(simulateGP)
library(TwoSampleMR)
#library(progress)


# Function allowing many simulations at a time in bc4
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
  
  # Adding the condition for replication: yes/no
  param <- cbind(param,repl = "y")
  
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
  dir.create("../results/sim", recursive=TRUE, showWarnings=FALSE)
  save(out, file=paste0("../results/sim/out", chunk, ".rdata"))
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
  
  # Adding the results to param
  param$b <- res$b
  param$se <- res$se
  param$nsnp_inc <- nrow(dat)
  param$pval <- res$pval
  param$mean_f <- mean(exp_d$fval, na.rm=TRUE)
  param$obs <- lm(y ~ x)$coef[2]
  param$what <- "all"
  
  # filter by most significant SNPs
  dat <- subset(dat, pval.exposure < 5e-8)
  
  # SNPs from the exposure replication dataset that gave significant results
  exp_d_s <- subset(exp_d, snp %in% dat$SNP)
  
  # SNPs from the exposure replication dataset that gave significant results
  exp_r_s <- subset(exp_r, snp %in% dat$SNP)
  
  # SNPs from the outcome dataset that correspond to the significant results
  out_s <- subset(out, snp %in% dat$SNP)
  
  # re-do the MR analysis for the significant SNPs, if there are any
  if (nrow(dat) > 0)
  {
    # specify the datasets for the MR analysis on the significant SNPs; the original dataset by default (no repl)
    exp = exp_d_s
    out = out_s
    # check if the replication option selected
    if (param$repl=="y")
    {
      exp = exp_r_s
    }
    
    dat <- data.frame(
      SNP = 1:nrow(exp),
      beta.exposure = exp$bhat,
      beta.outcome = out$bhat,
      se.exposure = exp$se,
      se.outcome = out$se,
      pval.exposure = exp$pval,
      pval.outcome = out$pval,
      fval.exposure = exp$fval,
      id.exposure=1,
      id.outcome=2,
      exposure=1,
      outcome=1,
      mr_keep=TRUE
    )
    res2 <- suppressMessages(mr(dat, method=c("mr_wald_ratio", "mr_ivw")))
    # Adding the results to p2
    p2$b <- res2$b
    p2$se <- res2$se
    p2$nsnp_inc <- nrow(dat)
    p2$pval <- res2$pval
    p2$mean_f <- mean(dat$fval.exposure, na.rm=TRUE)
    p2$obs <- lm(y ~ x)$coef[2]
  }
  p2$what <- "sig"
  
  # Saving the output from both sig and all
  out <- bind_rows(param, p2)
  return(out)
  
}

main()

