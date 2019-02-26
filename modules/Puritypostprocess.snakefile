#module: Purity PostProcessing steps by Facets in R
#import os
#from string import Template
_purityprocessing_threads=32


def purityprocessing_runsHelper(wildcards, iindex):
"""Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
  tmp = []
  r = config['runs'][wildcards.run]
  #print(r)

  # check that we have a valid pair
  if len(r) >=2:
  sample_name = r[iindex]
    tmp.append(input_template.format(sample=sample_name))
  else:
    #NOTE: I can't figure out a proper kill command so I'll do this
    tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
  return tmp

def getNormal_sample(wildcards):
   return purityprocessing_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
   return purityprocessing_runsHelper(wildcards, 1)

def purityprocessing_targets(wildcards):
   """Generates the targets for this module"""
   ls = []
   for run in config['runs']:
      ls.append("analysis/germlineVariants/%s/%s_purity_postprocessed_results.txt" % (run,run))
      ls.append("analysis/germlineVariants/%s/%s_results" % (run,run))
   return ls

rule purityprocessing_all:
 input:
   purityprocessing_targets

rule purityprocessing_filter:
 input:
    "analysis/purity/{run}/{run}_purity_results.txt"
 params:
    #no params
 output:
    "analysis/purity/{run}/{run}_purity_postprocessed_results.txt"
 benchmark:
    "benchmarks/puritycalls/{run}/{run}.postprocessed.purity.results.txt"
 shell:
    #"cidc_wes/modules/scripts/vcf_filterByReadDepth.py -v {input} -t {params.threshold} -f {params.field} -o {output}"
    """cat {input} | sed 's/chr//g'  > {output} """

rule purityplots_postprocessing:
 input:
   "analysis/purity/{run}/{run}_purity_postprocessed_results.txt"
 params:
   runname={run}
 output:
   output_dir="analysis/purity/{run}/{run}_results"
 benchmark:
    "benchmarks/puritycalls/{run}/{run}.purity.postprocessingplots.txt"
 shell:
    "Rscript --vanilla cidc_wes/modules/scripts/facets_plot.R {input} {output.output_dir} {params.runname}"
