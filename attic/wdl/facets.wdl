workflow Facets {
	call facets
}

task facets {
	File tumorBam
	File tumorBamIndex
	File normalBam
	File normalBamIndex
	File vcftar
	String pairName

	Int? diskSpaceGb = 80
	Int? memoryGb = 32
    String? docker = "gcr.io/cidc-biofx/facets:v2"

	command <<<

		# Remember starting directory b/c FireCloud doesn't let you put outputs of files w/ abs paths
		home_dir=$PWD

		cd /usr/gitc

		# Pileup mutations
        # REF: https://github.com/veseshan/facets/tree/master/inst/extcode
		./snp-pileup -q15 -Q20 ${vcftar} /tmp/${pairName} ${normalBam} ${tumorBam}

		# Remove chr in chromosome names
        cat /tmp/${pairName} | sed 's/chr//g' | gzip -c - > /usr/gitc/pileup/${pairName}.txt.gz

		# Gotta add a little vanilla
        echo '' > /usr/gitc/run_facets.R
        cat <<EOT >> /usr/gitc/run_facets.R

#!/usr/bin/env Rscript
setwd("/usr/gitc/pileup")
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

EOT
		Rscript --vanilla /usr/gitc/run_facets.R

		cd /usr/gitc/pileup

		# move to home directory for output
		mv /usr/gitc/pileup/Facets_output.txt $home_dir/${pairName}_Facets_output.txt
		mv /usr/gitc/pileup/Facets_iterations.txt $home_dir/Facets_iterations.txt
		mv /usr/gitc/pileup/${pairName}.txt.gz $home_dir/${pairName}_pileup.txt.gz
		mv /usr/gitc/pileup/${pairName}.txt.gz.cncf $home_dir/${pairName}_cncf.txt
		mv /usr/gitc/pileup/Rplots.pdf $home_dir/${pairName}_plot.pdf

		cd $home_dir
	>>>

	output {
		File facetsOutput = "${pairName}_Facets_output.txt"
        Array[Array[String]] facetsEstimation = read_tsv("${pairName}_Facets_output.txt")
		File facetsIterations = "Facets_iterations.txt"
		File pileup = "${pairName}_pileup.txt.gz"
		File cncf = "${pairName}_cncf.txt"
		File plot = "${pairName}_plot.pdf"
	}

	runtime {
		docker: "${docker}"
		slurm_docker: "${docker}"
		memory: "${memoryGb} GB"
		slurm_memory: "${memoryGb}G"
		cpu: "1"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}
