#!/usr/bin/env Rscript 
# @Author: Jingxin Fu
# @Email:  jingxinfu.tj@gmail.com
# @Create Date: 27 Dec 2018
# @License: MIT: https://opensource.org/licenses/MIT
# @R version: 3.5.1

######################################################################################################
# This Function 
# Call:
# tri_mutation_profile.R 
######################################################################################################

main <- function(){
  # load packages
  get_packages(requirements= c('optparse','MutationalPatterns','BSgenome','BSgenome.Hsapiens.UCSC.hg38'))

  # parsing parameters
  args <- parse_input()
  
}
tri_mutation_profile <- function(vcf_dir,out_dir) {
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
  vcf_files =  paste0(vcf_dir,dir(vcf_dir))
  sample_names = sapply(dir(vcf_dir),
                      function(x){
                        strsplit(x,split="\\.")[[1]][1]
                      })
  vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
  mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
}
parse_input <-function(){
  args <- commandArgs(trailingOnly=TRUE)
  # make option list and parse command line
  option_list <- list(
    make_option(c("-i", "--cellPath"), type="character", help="[Required]"),
    make_option(c("-c", "--cdfname"), type="character", default='a',help="Annotation CDF name [default %default]")
    make_option(c("-n", "--add_numbers"), action="store_true", default=FALSE,help="[default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list), args=args)
  # paramenter checking
  if(is.null(opts$cellPath)) stop('Location of .cel file must be provided. See script usage (--help)')

  return(opts)
}
get_packages <- function(requirements){
  suppressMessages(source("https://bioconductor.org/biocLite.R"))
  for(el in requirements){
    if (!is.element(el, installed.packages()[,1])){
      biocLite(el)
    }
    suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
  }
}

## call
if(!interactive()) {
   main()
}
