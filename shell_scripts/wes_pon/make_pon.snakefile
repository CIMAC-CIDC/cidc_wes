# Len Taing 2021 (TGBTG)
# Generate a panel of normal according to
# ref: https://support.sentieon.com/manual/TNscope_usage/tnscope/

from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

configfile: "make_pon.config.yaml"
def targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/align/%s/%s.sorted.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.bam.bai" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.dedup.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.dedup.bam.bai" % (sample,sample))
        ls.append("analysis/align/%s/%s.realigned.bam" % (sample,sample))
    	ls.append("analysis/align/%s/%s_prerecal_data.table" % (sample,sample))
    	ls.append("analysis/align/%s/%s_recalibrated.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_tnscope.output.vcf.gz" % (sample,sample))
    # ls.append("analysis/pon/%s.pon.vcf.gz" % config.get('center', 'cidc'))
    # ls.append("analysis/pon/%s.pon.vcf.gz.tbi" % config.get('center', 'cidc'))
    ls.append(GS.remote("%s/%s.pon.vcf.gz" % (config['destination'], config['center'])))
    ls.append(GS.remote("%s/%s.pon.vcf.gz.tbi" % (config['destination'], config['center'])))
    return ls

rule all:
    input:
        targets

def align_getFastq(wildcards):
    #ls = config["samples"][wildcards.sample]
    ls =[]
    for fastq in config['samples'][wildcards.sample]:
        ls.append(GS.remote(fastq))
    return ls

    return ls

def align_getBam(wildcards):
    bam = config["samples"][wildcards.sample][0] #returns only the first elm
    return bam

#BASED on xindong's aggregate_align_input in CIDC_Chips
def aggregate_align_input(wildcards):
    # handle .bam files separately from .fastq files
    #check only the first file
    sample_first_file = config["samples"][wildcards.sample][0]
    if sample_first_file.endswith(".bam"):
        return ["analysis/align/{sample}/{sample}.sorted.fromBam.bam",
                "analysis/align/{sample}/{sample}.sorted.fromBam.bam.bai"]
    else:
        return ["analysis/align/{sample}/{sample}.sorted.fromFastq.bam",
                "analysis/align/{sample}/{sample}.sorted.fromFastq.bam.bai"]

rule aggregate_input:
    input:
        aggregate_align_input
    params:
        bam = lambda wildcards,input: input[0],
        bai = lambda wildcards,input: input[1],
    output:
       bam="analysis/align/{sample}/{sample}.sorted.bam",
       bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
    shell:
        "mv {params.bam} {output.bam} && mv {params.bai} {output.bai}"

rule align_from_bam:
    input:
        align_getBam
    output:
        bam="analysis/align/{sample}/{sample}.sorted.fromBam.bam",
        bai="analysis/align/{sample}/{sample}.sorted.fromBam.bam.bai"
    threads: 32 #_bwa_threads
    priority: 100
    params:
        sentieon_path=config['sentieon_path'],
        bwa_index=config['bwa_index'],
        #DON'T write a mini-program withi a program-
        #awk cmd to add sample names to RGs!!
        awk_cmd=lambda wildcards: "awk -v OFS=\'\\t\' \'{ split($2,a,\":\"); read_id=a[2]; $2=\"ID:%s.\" read_id; gsub(/SM:.+\\t/,\"SM:%s\\t\"); print $0}\'" % (wildcards.sample, wildcards.sample),
        #NEVER do it twice!- gawk cmd to inject sample name into each read!!!
        gawk_cmd=lambda wildcards: "gawk -v OFS=\'\\t\' \'{rg=match($0,/RG:Z:(\S+)/,a); read_id=a[1]; if (rg) {sub(/RG:Z:\S+/, \"RG:Z:%s.\" read_id, $0); print $0} else { print $0 }}\'" % wildcards.sample,
    benchmark: "benchmarks/align/{sample}/{sample}.align_from_bam.txt"
    shell:
        """samtools view -H {input} | grep \"^@RG\" | {params.awk_cmd} > {wildcards.sample}.header && samtools collate --output-fmt SAM -@ {threads} -Of {input} | {params.gawk_cmd} | samtools view -@ {threads} -b - | samtools fastq -@ {threads} -t -s /dev/null -0 /dev/null - | ({params.sentieon_path}/sentieon bwa mem -t {threads} -M -K 10000000 -p -C -H {wildcards.sample}.header {params.bwa_index} - || echo -n 'error' ) | {params.sentieon_path}/sentieon util sort -t {threads} -o {output.bam} --sam2bam -; rm {wildcards.sample}.header"""

rule align_from_fastq:
    input:
        align_getFastq
    output:
        bam="analysis/align/{sample}/{sample}.sorted.fromFastq.bam",
        bai="analysis/align/{sample}/{sample}.sorted.fromFastq.bam.bai"
    params:
        sentieon_path=config['sentieon_path'],
        bwa_index=config['bwa_index'],
        read_group= lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample),        
        input_bases="10000000",
        #need to adjust threads for the other process
        tthreads=lambda wildcards, input, output, threads, resources: threads-1
    threads: 32 #_bwa_threads
    priority: 100
    message: "ALIGN: Running sentieon BWA mem for alignment"
    log: "analysis/logs/align/{sample}/align.sentieon_bwa.{sample}.log"
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.align_from_fastq.txt"
    shell:
        """({params.sentieon_path}/sentieon bwa mem -M -R \"{params.read_group}\" -t {params.tthreads} -K {params.input_bases} {params.bwa_index} {input} || echo -n 'error' ) | {params.sentieon_path}/sentieon util sort -r {params.bwa_index} -o {output.bam} --sam2bam -i -"""

rule scoreSample:
    "Calls sentieon driver  --fun score_info on the sample"
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
    output:
        score="analysis/align/{sample}/{sample}.sorted.score.txt",
        idx="analysis/align/{sample}/{sample}.sorted.score.txt.idx",
    message: "ALIGN: score sample"
    log: "analysis/logs/align/{sample}/align.scoreSample.{sample}.log"
    threads: 8 #_align_threads
    params:
        index1=config['sentieon_path'],
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.scoreSample.txt"
    shell:
        """{params.index1}/sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output.score}"""

rule dedupSortedUniqueBam:
    """Dedup sorted unique bams using sentieon
     output {sample}_unique.sorted.dedup.bam"""
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
        score="analysis/align/{sample}/{sample}.sorted.score.txt"
    output:
        bamm="analysis/align/{sample}/{sample}.sorted.dedup.bam",
        baii="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai",
        met="analysis/align/{sample}/{sample}.sorted.dedup.metric.txt",
    message: "ALIGN: dedup sorted unique bam file"
    log: "analysis/logs/align/{sample}/align.dedupSortedUniqueBam.{sample}.log"
    threads: 32 #_align_threads
    priority: 75
    params:
        index1=config['sentieon_path'],
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.dedupSortedUniqueBam.txt"
    shell:
        """{params.index1}/sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {input.score} --metrics {output.met} {output.bamm}"""

#==============================================================================

rule Indel_realigner_sentieon:
    """indel realigner for uniquely  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.dedup.bam",
         bai="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai",
    output:
         realignbam="analysis/align/{sample}/{sample}.realigned.bam",
         realignbai="analysis/align/{sample}/{sample}.realigned.bam.bai"
    message:
         "INDEL REALIGNER: indel realigner for  mapped reads"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        #dbsnp=config['dbsnp'], #not used!
        mills=config['Mills_indels'],
        g1000=config['G1000_indels'],
    group: "recalibration"
    threads: 8 #_realigner_threads
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Indel_realigner_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.bam} --algo Realigner -k {params.mills} -k {params.g1000} {output.realignbam}"""

rule Base_recalibration_precal_sentieon:
    """base recalibration for realigned files"""
    input:
        realignbam="analysis/align/{sample}/{sample}.realigned.bam",
    output:
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table",
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",
        recalibratedbai="analysis/align/{sample}/{sample}_recalibrated.bam.bai"
    message:
        " PRE BASE RECALIBRATION: base recalibration for  realigned files"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: 4 #_realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_precal_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.realignbam} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000}  {output.prerecaltable} --algo ReadWriter {output.recalibratedbam}"""
#==============================================================================

rule somatic_tnscope:
    input:
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",
        recalibratedbai="analysis/align/{sample}/{sample}_recalibrated.bam.bai"
    output:
        tnscopevcf="analysis/align/{sample}/{sample}_tnscope.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        #pon= config['pon'],
        #dbsnp= config['dbsnp'],
        normal = lambda wildcards: wildcards.sample,
    threads: 18 #_somatic_threads
    priority: 50
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{sample}/{sample}.somatic_calling_TNscope.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input.recalibratedbam} --algo TNscope --tumor_sample {params.normal} {output.tnscopevcf}"""

rule mergeVCFs:
    input:
        expand("analysis/align/{sample}/{sample}_tnscope.output.vcf.gz", sample=config['samples'])
    output:
        "analysis/pon/{center}.pon.vcf"
    benchmark:
        "benchmarks/mergeVCFs.{center}.txt"
    shell:
        "bcftools merge -m all -f PASS,. --force-samples {input} | "
        "bcftools plugin fill-AN-AC | "
        "bcftools filter -i 'SUM(AC)>1' > {output}"

rule gzipAndTabix:
    input:
        "analysis/pon/{center}.pon.vcf"
    output:
        gz=GS.remote("%s/{center}.pon.vcf.gz" % (config['destination'])),
	tbi=GS.remote("%s/{center}.pon.vcf.gz.tbi" % (config['destination']))
        # gz="analysis/pon/{center}.pon.vcf.gz",
        # tbi="analysis/pon/{center}.pon.vcf.gz.tbi",
    benchmark:
        "benchmarks/gzipAndTabix.{center}.txt"
    shell:
        "bgzip -c {input} > {output.gz} && tabix -p vcf {output.gz}"
