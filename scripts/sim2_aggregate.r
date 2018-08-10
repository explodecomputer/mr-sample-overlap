library(dplyr)

fn <- list.files("../results/sim2", full=TRUE)
l <- list()
for(i in 1:length(fn))
{
	message(i)
	load(fn[i])
	l[[i]] <- out
}

param <- bind_rows(l)
save(param, file="../results/sim2.rdata")



