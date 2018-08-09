library(simulateGP)
library(TwoSampleMR)
library(progress)

main <- function()
{
	param <- expand.grid(
		nsim = 1:100,
		nid=c(10000),
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
	set.seed(chunk)

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
	dir.create("../results/scenario1", recursive=TRUE, showWarnings=FALSE)
	save(out, file=paste0("../results/scenario1/out", chunk, ".rdata"))
}

get_r_vals <- function(f, n)
{
	sqrt(f / (n-2+f))
}


runsim <- function(param)
{
	g1 <- make_geno(param$nid, 1, 0.5)
	g2 <- make_geno(param$nid, 1, 0.5)
	u1 <- rnorm(param$nid)
	u2 <- rnorm(param$nid)
	gx <- get_r_vals(param$f, param$nid)
	x1 <- make_phen(c(gx, param$ux), cbind(g1, u1))
	x2 <- make_phen(c(gx, param$ux), cbind(g2, u2))
	y1 <- make_phen(c(param$xy, param$uy), cbind(x1, u1))
	y2 <- make_phen(c(param$xy, param$uy), cbind(x2, u2))
	overlap_index1 <- c(1:floor(param$nid * param$overlap))
	overlap_index2 <- c(1:param$nid)[! c(1:param$nid) %in% overlap_index1]
	go <- rbind(g1[overlap_index1,, drop=FALSE], g2[overlap_index2,, drop=FALSE])
	xo <- c(x1[overlap_index1], x2[overlap_index2])
	yo <- c(y1[overlap_index1], y2[overlap_index2])


	obs11 <- fast_assoc(x1, y1)
	obs22 <- fast_assoc(x2, y2)
	obso <- fast_assoc(x1, yo)
	param$obs11_beta <- obs11$bhat
	param$obs11_se <- obs11$se
	param$obs11_pval <- obs11$pval
	param$obs22_beta <- obs22$bhat
	param$obs22_se <- obs22$se
	param$obs22_pval <- obs22$pval
	param$obso_beta <- obso$bhat
	param$obso_se <- obso$se
	param$obso_pval <- obso$pval

	dat11 <- get_effs(x1, y1, g1)
	dat22 <- get_effs(x2, y2, g2)
	dato <- get_effs(xo, yo, go)
	names.exposure <- c(1, grep("exposure", names(dat11)), 14)
	names.outcome <- grep("outcome", names(dat11))
	dat12 <- cbind(dat11[,names.exposure], dat22[,names.outcome])
	dat21 <- cbind(dat22[,names.exposure], dat11[,names.outcome])
	dato <- cbind(dat11[,names.exposure], dato[,names.outcome])
	mr11 <- mr(dat11)
	mr22 <- mr(dat22)
	mr12 <- mr(dat12)
	mr21 <- mr(dat21)
	mro <- mr(dato)

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
	param$mro_beta <- mro$b
	param$mro_se <- mro$se
	param$mro_pval <- mro$pval

	param$fval11 <- fast_assoc(g1, x1)$fval
	param$fval22 <- fast_assoc(g2, x2)$fval
	param$fvalo <- fast_assoc(go, xo)$fval
	return(param)
}

main()
