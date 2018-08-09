library(simulateGP)
library(TwoSampleMR)
library(progress)

main <- function()
{
	param <- expand.grid(
		nsim = 1:100,
		nid1=c(10000),
		nid2=c(10000),
		overlap=seq(0, 1, by=0.2),
		f = c(2,4,6,8,10,12,15, 20, 30, 50, 100),
		xy = c(0, 0.01, 0.1, -0.01, -0.1),
		ux = c(0, -0.01, -0.1, 0.01, 0.1),
		uy = c(0, -0.01, -0.1, 0.01, 0.1)
	)
	# param <- sc1[555500,]

	args <- commandArgs(T)
	splitsize <- as.numeric(args[1])
	chunk <- as.numeric(args[2])
	first <- (chunk-1) * splitsize + 1
	last <- min(chunk * splitsize, nrow(param))
	param <- param[first:last,]
	message("Running ", nrow(param), " simulations")
	message(first, " to ", last)

	out <- list()
	pb <- progress_bar$new(total=nrow(param))
	for(i in 1:nrow(param))
	{
		pb$tick()
		out[[i]] <- suppressMessages(runsim(param[i,]))
	}
	out <- bind_rows(out)
	dir.create("../results/scenario1", recursive=TRUE)
	save(out, file=paste0("../results/scenario1/out", chunk, ".rdata"))
}

get_r_vals <- function(f, n)
{
	sqrt(f / (n-2+f))
}


runsim <- function(param)
{
	g1 <- make_geno(param$nid1, 1, 0.5)
	g2 <- make_geno(param$nid2, 1, 0.5)
	u1 <- rnorm(param$nid1)
	u2 <- rnorm(param$nid2)
	gx <- get_r_vals(param$f, param$nid1)
	x1 <- make_phen(c(gx, param$ux), cbind(g1, u1))
	x2 <- make_phen(c(gx, param$ux), cbind(g2, u2))
	y1 <- make_phen(c(param$xy, param$uy), cbind(x1, u1))
	y2 <- make_phen(c(param$xy, param$uy), cbind(x2, u2))

	obs11 <- fast_assoc(x1, y1)
	obs22 <- fast_assoc(x2, y2)
	param$obs11_beta <- obs11$bhat
	param$obs11_se <- obs11$se
	param$obs11_pval <- obs11$pval
	param$obs22_beta <- obs22$bhat
	param$obs22_se <- obs22$se
	param$obs22_pval <- obs22$pval

	dat11 <- get_effs(x1, y1, g1)
	dat22 <- get_effs(x2, y2, g2)
	names.exposure <- c(1, grep("exposure", names(dat11)), 14)
	names.outcome <- grep("outcome", names(dat11))
	dat12 <- cbind(dat11[,names.exposure], dat22[,names.outcome])
	dat21 <- cbind(dat22[,names.exposure], dat11[,names.outcome])
	mr11 <- mr(dat11)
	mr22 <- mr(dat22)
	mr12 <- mr(dat12)
	mr21 <- mr(dat21)

	param$mr11_beta <- mr11$b
	param$mr11_se <- mr11$se
	param$mr11_pval <- mr11$pval
	param$mr22_beta <- mr22$b
	param$mr22_se <- mr22$se
	param$mr22_pval <- mr22$pval
	param$mr12_beta <- mr12$b
	param$mr12_se <- mr12$se
	param$mr12_pval <- mr12$pval
	param$mr21_beta <- mr21$b
	param$mr21_se <- mr21$se
	param$mr21_pval <- mr21$pval

	param$fval11 <- fast_assoc(g1, x1)$fval
	param$fval22 <- fast_assoc(g2, x2)$fval
	return(param)
}

main()
