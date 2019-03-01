<<<<<<< HEAD
#!/usr/bin/env Rscript  
#setRepositories(graphics = FALSE, ind = 1:6)
#install.packages("sequenza")
#install.packages("devtools")
library(devtools)
install_github("aroneklund/copynumber")
library(sequenza)

sequenza_results<-function(arg_in,arg_out,arg_name){
    test <- sequenza.extract(arg_in,assembly="hg38")
=======
#!/usr/bin/env Rscript
#setRepositories(graphics = FALSE, ind = 1:6)
#install.packages("sequenza")

library(sequenza)

sequenza_results<-function(arg_in,arg_out,arg_name){
    test <- sequenza.extract(arg_in)
>>>>>>> 33a19477feb624e0f31c2c100e065da78025f190
    CP <- sequenza.fit(test)
    sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = arg_name, out.dir=arg_out)
    mut.tab <- read.table(paste0(arg_out,"/",arg_name,"_mutations.txt"), header=TRUE)
    seg.res <- read.table(paste0(arg_out,"/",arg_name,"_segments.txt"),header=TRUE)
    ## fix chr list. Remark: do not use '$'
    chrs = names(test[['mutations']])
    mut.tab[['chromosome']] = factor(as.character(mut.tab[['chromosome']]), levels=chrs)
    seg.res[['chromosome']] = factor(as.character(seg.res[['chromosome']]), levels=chrs)

    tsv <- sequenza:::sequenza2PyClone(mut.tab, seg.res, "arg_name", norm.cn = 2)
    tsv <- tsv[tsv[,'major_cn'] != 0 ,] # Required for Pyclone
    write.table(tsv,paste0(arg_out,"/",arg_name,'.pyclone.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
}

args <- commandArgs( trailingOnly = TRUE )
arg_in= args[1]
arg_out = args[2]
arg_name = args[3]

sequenza_results(arg_in,arg_out,arg_name)
<<<<<<< HEAD
=======





>>>>>>> 33a19477feb624e0f31c2c100e065da78025f190
