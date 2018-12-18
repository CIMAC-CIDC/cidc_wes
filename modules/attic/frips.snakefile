#MODULE: FRiPs- Calculating the fraction of reads under peak (score/percentage)

#helper fns to get samples associated for each run
#copied directly from peaks.snakefile
#TODO: centralize these helper fns!

#PARAMETERS:
_logfile="analysis/logs/frips.log"
_macs_fdr="0.01"
_macs_keepdup="1"
_macs_extsize="146"
_macs_species="hs"
_samtools_threads=4

#NOTE: using the _refs from chips.snakefile
def frips_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/frag/%s/%s_fragDist.png" % (sample,sample))
        ls.append("analysis/frips/%s/%s_pbc.txt" % (sample,sample))
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            ls.append("analysis/frips/%s/%s_frip.txt" % (runRep,runRep))
    ls.append("analysis/frips/pbc.csv")
    ls.append("analysis/frips/nonChrM_stats.csv")
    ls.append("analysis/frips/frips.csv")
    return ls

#NOTE: this takes the Treatment of the FIRST replicate!!!--ignores the rest!
def frip_getTreatBam(wildcards):
    """RETURNS the associated 4M_nonChrM.bam for the run's treatment sample"""
    r = config['runs'][wildcards.run]
    #GET first treatement sample
    first = r[0]
    #print(first)
    ret = "analysis/align/%s/%s_4M_nonChrM.bam" % (first, first)
    return ret

rule frips_all:
    input:
        frips_targets

#SECTION: generate 4M sub-sampled bams
#It's a bit redundant, but I tried my best to reduce the redundancy
rule create_nonChrM:
    """RULE that will generate the base file for both 
    4M_nonChrM.bam AND 4M_unique_nonChrM.bam so we save some redundancy here"""
    input:
        "analysis/align/{sample}/{sample}.sorted.bam"
    params:
        #hack to get the regex in to filter out chrM, random, chrUn
        regex="\'/chrM/d;/random/d;/chrUn/d\'",
        #msg just for message below
        msg= lambda wildcards: wildcards.sample
    message: "FRiPs: creating the nonChrM SAM file {params.msg}"
    log:_logfile
    threads: _samtools_threads
    output:
        temp('analysis/align/{sample}/{sample}_nonChrM.sam')
    shell:
        "samtools view -@ {threads} -h {input} | sed -e {params.regex} > {output} 2>>{log}"

rule sample_4M_from_nonChrM:
    """Sample 4M reads from nonChrM SAM file (from create nonChrM)
    ref: https://sourceforge.net/p/samtools/mailman/message/29011091/ 
    see '-s 21.5'
    """
    input:
        'analysis/align/{sample}/{sample}_nonChrM.sam'
    params:
        n="4000000"
    message: "FRiPs: sample- 4M from non-chrM reads"
    log:_logfile
    output:
        temp('analysis/align/{sample}/{sample}_4M_nonChrM.bam')
    shell:
        """
        chips/modules/scripts/frips_sample.sh -n {params.n} -i {input} -o {output} 2>>{log}
        """

rule create_unique_nonChrM:
    """Generate _unique_nonChrM.bam by
    Filter out non-uniquely mapped reads from _nonChrM.sam
    """
    input:
        "analysis/align/{sample}/{sample}_nonChrM.sam"
    params:
        #msg just for message below
        msg= lambda wildcards: wildcards.sample
    message: "FRiPs: create uniquely mapped non-chrM reads {params.msg}"
    log:_logfile
    threads: _samtools_threads
    output:
        temp('analysis/align/{sample}/{sample}_unique_nonChrM.bam')
    shell:
        "samtools view -@ {threads} -b -h -F 4 {input} > {output} 2>>{log}"

rule sample_4M_from_uniqueNonChrM:
    """Sample 4M reads from uniqueNonChrM reads
    ref: https://sourceforge.net/p/samtools/mailman/message/29011091/ 
    see '-s 21.5'
    """
    input:
        'analysis/align/{sample}/{sample}_unique_nonChrM.bam'
    params:
        n="4000000",
        #msg just for message below
        msg= lambda wildcards: wildcards.sample
    message: "FRiPs: sample- 4M from uniquely mapped non-chrM reads {params.msg}"
    log:_logfile
    output:
        temp('analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam')
    shell:
        """
        chips/modules/scripts/frips_sample.sh -n {params.n} -i {input} -o {output} 2>>{log}
        """

rule frip_calculate:
    """Calculate the frip score"""
    #TODO: if there are more than 1 treatment, merge them??!
    input:
        treat=frip_getTreatBam,
        bed="analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak",
    output:
        "analysis/frips/{run}.{rep}/{run}.{rep}_frip.txt"
    params:
        pval="1E-9"
    message: "FRiPs: calculate frips"
    log:_logfile
    shell:
        "chips/modules/scripts/frips_calculate.sh -a {input.treat} -b {input.bed} -p {params.pval} > {output} 2>>{log}"

rule frip_pbc:
    """Generate the PBC histogram for each normalized sample, which will be 
    used to calculate N1, Nd, and PBC (for the report)
    """
    input:
        "analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam"
    output:
        "analysis/frips/{sample}/{sample}_pbc.txt"
    params:
        #msg just for message below
        msg= lambda wildcards: wildcards.sample
    message: "FRiP: generate PBC histogram for each sample/bam {params.msg}"
    log: _logfile
    shell:
        "chips/modules/scripts/frips_pbc.sh -i {input} -o {output} 2>> {log}"

rule collect_pbc:
    """Collect and parse out the PBC for the ALL of the samples"""
    input:
        expand("analysis/frips/{sample}/{sample}_pbc.txt", sample=sorted(config["samples"]))
    output:
        "analysis/frips/pbc.csv"
    message: "ALIGN: collect and parse ALL pbc stats"
    log: _logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frips_collectPBC.py -f {files} > {output} 2>>{log}")

rule nonChrM_stats:
    """Get the nonChrM mapping stats for each aligment run"""
    input:
        #NOTE: uniq_bam is generated in align_common module-
        #HACK- taking advantage that this module is loaded AFTER align_common
        uniq_bam="analysis/align/{sample}/{sample}_unique.bam",
        nonChrM_bam="analysis/align/{sample}/{sample}_unique_nonChrM.bam"
    output:
        #NOTE: CHIPS-RESUME- IT's better to keep this file around, b/c
        #it's inputs are both temp files that are hard to generate
        "analysis/align/{sample}/{sample}_nonChrM_stat.txt"
    message: "ALIGN: get nonChrM mapping stats for each bam"
    log: _logfile
    params:
        sam_th = _samtools_threads / 2
    threads: _samtools_threads
    shell:
        #FLAGSTATS is the top of the file, and we append the uniquely mapped
        #reads to the end of the file
        "samtools view -@ {params.sam_th} -c {input.uniq_bam} > {output} 2>>{log}"
        " && samtools view -@ {params.sam_th} -c {input.nonChrM_bam} >> {output} 2>> {log}"

rule collect_nonChrM_stats:
    """Aggregate all nonChrM stats for ALL of the samples"""
    input:
        expand("analysis/align/{sample}/{sample}_nonChrM_stat.txt", sample=sorted(config["samples"]))
    output:
        "analysis/frips/nonChrM_stats.csv"
    message: "FRiPs: collect and parse ALL nonChrM stats"
    log: _logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frips_collectNonChrM.py -f {files} > {output} 2>>{log}")

rule get_SampleFragLength:
    """Dump all of the sample's fragment lengths into 
    frag/{sample}/{sample}_frags.txt, so we can generate the distribution plot in make_FragPlot
    """
    input:
        "analysis/align/{sample}/{sample}_unique.sorted.bam"
    params:
        awk_cmd = """awk ' $1 <= 1000 && $1 > 0 '"""
    message: "FRAG: get fragment sizes"
    log:_logfile
    threads: _samtools_threads
    output:
        temp("analysis/frag/{sample}/{sample}_frags.txt")
    shell:
        #GRAB out the 9th column, ensuring it's in 1-1000
        "samtools view -@ {threads} {input} | cut -f 9 | {params.awk_cmd} > {output} 2>>{log}"
 
rule make_FragPlot:
    """plot the fragment distribution:
    generate the R plot by running frag_plotFragDist.R on _frags.txt
    """
    input:
        "analysis/frag/{sample}/{sample}_frags.txt"
    params:
        name= lambda wildcards: wildcards.sample
    message: "FRAG: plot fragment size distribution plot"
    log:_logfile
    output:
        "analysis/frag/{sample}/{sample}_fragDist.png"
    shell:
        #RUN the R script to get the plot
        "chips/modules/scripts/frag_plotFragDist.R {input} {output} {params.name} 2>>{log}"
    
rule getFripStats:
    """Collect the frips statistics from analysis/frips/{run}/{run}_frip.txt"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput("analysis/frips/$runRep/$runRep_frip.txt")
    output:
        "analysis/frips/frips.csv"
    message: "FRiPs: collecting frips stats for each run"
    log:_logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/frips_getFrips.py -f {files} -o {output} 2>>{log}")
