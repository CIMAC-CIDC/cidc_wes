#!/usr/bin/env Rscript                                                                                                                                                                       #!/usr/bin/env Rscript
library(pctGCdata)
library(facets)

facets_plots<-function(arg_in,arg_out, arg_name) {
  rcmat = readSnpMatrix(arg_in, perl.pileup = F)
  #print(rcmat)
  #stop()
  filenames <- c()
  modes <- c()
  puritys <- c()
  ploidys <- c()
  dipLogRs <- c()
  purity.vals <- c()
  ploidy.vals <- c()
  dipLogR.vals <- c()
  cncf.vals <- list()

  getmode <- function(v) {
     uniqv <- unique(v)
     uniqv[which.max(tabulate(match(v, uniqv)))]
   }


  for(i in seq(10)) {
      xx = preProcSample(rcmat)
      oo=procSample(xx, cval=150)
      fit=emcncf(oo, maxiter = 1000)
      purity.vals[i] <- signif(fit[['purity']], 3)
      ploidy.vals[i] <- fit[['ploidy']]
      dipLogR.vals[i] <- oo[['dipLogR']]
      cncf.vals[[i]] <- fit[['cncf']]
      plotSample(x=oo, emfit=fit)
   }

   mode <- getmode(purity.vals)
   index <- which(purity.vals == mode)[1]
   purity <- purity.vals[index]
   ploidy <- ploidy.vals[index]
   dipLogR <- dipLogR.vals[index]
   write.table(data.frame(cncf.vals[[index]]), paste0(arg_out,arg_name,'.cncf'), quote=F, sep="\t", row.names=F)

   filenames <- c(filenames, arg_in)
   modes <- c(modes, mode)
   puritys <- c(puritys, purity)
   ploidys <- c(ploidys, ploidy)
   dipLogRs <- c(dipLogRs, dipLogR)
#}

  df <- data.frame(name=filenames, mode=modes, purity=puritys, ploidy=ploidys, dipLogR=dipLogRs)
  df.iter <- data.frame(name=filenames, purity=purity.vals, ploidy=ploidy.vals, dipLogR=dipLogR.vals)
  write.table(df,paste0(arg_out,arg_name,'.optimalpurityvalue.txt'), quote=F, sep="\t", row.names=F)
  write.table(df.iter, paste0(arg_out,arg_name,'.iterpurityvalues.txt') , quote=F, sep="\t", row.names=F)
}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
arg_out = args[2]
arg_name = args[3]

facets_plots(arg_in,arg_out,arg_name)


