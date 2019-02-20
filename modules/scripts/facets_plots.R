#!/usr/bin/env Rscript
setwd("/usr/gitc/pileup") ##set this to  the actual  pileup path-ask len 
library(pctGCdata)
library(facets)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

filenames <- c()
modes <- c()
puritys <- c()
ploidys <- c()
dipLogRs <- c()

for(file in list.files(pattern='.gz$')) {

  purity.vals <- c()
  ploidy.vals <- c()
  dipLogR.vals <- c()
  cncf.vals <- list()
  rcmat = readSnpMatrix(file, perl.pileup = F)

  for(i in seq(10)) {
    xx = preProcSample(rcmat)
    oo=procSample(xx, cval=150)
    fit=emcncf(oo, maxiter = 1000)
    purity.vals[i] <- signif(fit[['purity']], 3)
    ploidy.vals[i] <- fit[['ploidy']]
    dipLogR.vals[i] <- oo[['dipLogR']]
    cncf.vals[[i]] <- fit[['cncf']]
    plotSample(x=oo, emfit=fit) # plot
  }

  print(file)
  print(purity.vals)

  mode <- getmode(purity.vals)
  index <- which(purity.vals == mode)[1]
  purity <- purity.vals[index]
  ploidy <- ploidy.vals[index]
  dipLogR <- dipLogR.vals[index]
  write.table(data.frame(cncf.vals[[index]]), paste0(file,'.cncf'), quote=F, sep="\t", row.names=F)

  filenames <- c(filenames, file)
  modes <- c(modes, mode)
  puritys <- c(puritys, purity)
  ploidys <- c(ploidys, ploidy)
  dipLogRs <- c(dipLogRs, dipLogR)
}

df <- data.frame(name=filenames, mode=modes, purity=puritys, ploidy=ploidys, dipLogR=dipLogRs)
df.iter <- data.frame(name=filenames, purity=purity.vals, ploidy=ploidy.vals, dipLogR=dipLogR.vals)
write.table(df, "Facets_output.txt", quote=F, sep="\t", row.names=F)
write.table(df.iter, "Facets_iterations.txt", quote=F, sep="\t", row.names=F)
