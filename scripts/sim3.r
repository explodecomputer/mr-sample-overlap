library(simulateGP)
library(TwoSampleMR)
#library(progress)


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
	dir.create("../results/sim3", recursive=TRUE, showWarnings=FALSE)
	save(out, file=paste0("../results/sim3/out", chunk, ".rdata"))
}


get_r_vals <- function(f, n)
{
	sqrt(f / (n-2+f))
}


runsim <- function(param)
{
	p2 <- param
	g <- make_geno(param$nid, param$nsnp, param$maf)
	u <- rnorm(param$nid)
	x <- cbind(g,u) %*% c(rep(param$gx, param$nsnp), param$ux) + rnorm(param$nid)
	y <- param$xy * x + u * param$uy + rnorm(param$nid)

	exp <- gwas(x[1:(param$nid/2)], g[1:(param$nid/2),])
	out <- gwas(y[1:(param$nid/2+param$offset)], g[1:(param$nid/2+param$offset),])
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
	res <- suppressMessages(mr(dat, method=c("mr_wald_ratio", "mr_ivw")))

	param$b <- res$b
	param$se <- res$se
	param$nsnp_inc <- nrow(dat)
	param$pval <- res$pval
	param$mean_f <- mean(exp$fval, na.rm=TRUE)
	param$obs <- lm(y ~ x)$coef[2]
	param$what <- "all"

	dat <- subset(dat, pval.exposure < 5e-8)
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


main()

