# WES- Whole Exome-seq analysis pipeline in snakemake

# Introduction to WES

# Table of Contents

### Panel of normal for CNV 

We use the normal samples from Pilot 1(Broad and MDACC) and [Miao D, et al., Science, 2018](http://science.sciencemag.org/content/359/6377/801.abstract) . As Sentieon only accpets one target BED file to create the panel of normal, we intercept all BED files (only Broad and MDACC available) to generate the `target_bed`. The `target_padding` was 0 according to information from Broad. Following are cmd lines to generate the PoN:

```shell
#1. Create a panel of normal
echo "1.0 ----------Generate overlapped region among input BED files----------"
bedtools intersect -a $MDA_TARGET_BED -b $Broad_TARGET_BED | \
awk '$1 ~ "[XY0-9]$" {print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3}' > $OUT_PON_Dir/overlap.bed

echo "1.1 ----------Create a panel of normal for CNV----------"
$release_dir/bin/sentieon driver -t $NUMBER_THREADS -r $REFERENCE \
$( ls $Broad_normal_bams | while read s;do echo "-i ${s}";done ) \
$( ls $MDA_normal_bams | while read s;do echo "-i ${s}";done ) \
$( ls $Miao_normal_bams | while read s;do echo "-i ${s}";done ) \
--algo CNV \
--target $OUT_PON_Dir/overlap.bed \
--target_padding 0 \
--create_pon $OUT_PON_Dir/MDA_Broad_Miao_pad0
```



# Installing WES
You will only need to install WES once, either for your own use, or if you are a system administrator, for the entire system (see **Appendix C**).  In other words, you will only need to perform the steps described in this section only once.  
NOTE: this section ends with **Using WES** (below)

### Required software
We assume that the following tools are already installed on your system and that you have some basic familiarity in using them:
`git`
`wget`
### Installing Miniconda
WES uses the [Conda](https://conda.io/docs/intro.html) packaging system to manage and install all of its required software packages.
To install miniconda:
1.  download the Miniconda installer: 
    `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`
2.  run the installer:
    `bash Miniconda3-latest-Linux-x86_64.sh`
### Installing the WES conda environment
Conda environments are briefly explained [here](https://conda.io/docs/using/envs.html).  Briefly, if you are familiar with [Python Virtual Environments](http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/) or [Docker Containers](https://www.docker.com/what-container) then Conda environments should be a familiar concept.  

If you are **not familiar** with these concepts, then a conda environment is simply a **self-contained package space that is composed of various packages.**  So for example, a **bioinformatics** conda space may include packages such as **R**, **samtools**, **bedtools**, etc.

WES is dependent on the *wes* conda environment.
0. **clone the wes source code**:
    `git clone git@bitbucket.org:plumbers/cidc_wes
    ** NOTE: this command will create a directory called 'cidc_wes'.  After the next five steps, this directory can be safely deleted as we will explain how to *Setup a WES Project* below. **

1. **installing wes (conda env)**:
    After cloning the git repository, create the wes environment by doing this:
    - `cd cidc_wes`
    - `conda env create -f environment.yml`

### Downloading the Chips static reference files-
You will have to try to find the hg38 bwa index and fasta file

# Using WES
### Anatomy of a WES project
All work in WES is done in a **PROJECT** directory, which is simply a directory to contain a single WES analysis run.  **PROJECT** directories can be named anything (and they usually start with a simple mkdir command, e.g. mkdir wes_for_paper),  but what is CRITICAL about a **PROJECT** directory is that you fill them with the following core components:
(We first lay out the directory structure and explain each element below)
> PROJECT/
> cidc_wes/
> data/  - *optional*
> config.yaml
> metasheet.csv
> ref_dir

The 'cidc_wes' directory contains all of the WES source code.  We'll explain how to download that directory below.  The 'data' directory is an optional directory that contains all of your raw data.  It is optional because those paths __may__ be fully established in the config.yaml, __however__ it is best practice to gather your raw data within 'data' using [symbolic links](https://www.cyberciti.biz/faq/creating-soft-link-or-symbolic-link/).

The *config.yaml* and *metasheet.csv* are configurations for your WES run (also explained below).

The cidc_wes/ref.yaml file (not shown or discussed) is explained in **Appendix E**.

After a successful **WES** run, another 'analysis' folder is generated which contains all of the resulting output files.

### Setting up a WES project
0. **creating the project directory**
    As explained above, the **PROJECT** directory is simply a directory to contain an entire WES run.  **It can be named anything, but for this section, we'll simply call it 'PROJECT'**
    `mkdir PROJECT`
    `cd PROJECT`
1. **link data files**
    As explained above, creating a data directory is a place to gather all of your **raw data files (.fastq, .fastq.gz, .bam)**.  It is optional, but **highly recommended**.
    `mkdir data`
    And in 'data', copy over or make symbolic links to your raw data files
2. **clone cidc_wes**
    In your PROJECT directory:
    `git clone git@bitbucket.org:plumbers/cidc_wes`
3. **creating config.yaml and metasheet.csv**
    a. **copy cidc_wes/config.yaml and cidc_wes/metasheet.csv into the PROJECT dir:**
    In the PROJECT directory:
    `cp cidc_wes/config.yaml .`
    `cp cidc_wes/metasheet.csv .`

    b. **setup config.yaml**:
        The config.yaml is where you define WES run parameters and the ChIP-seq samples for analysis.
        

    1. **Set the assembly**: typically hg38 or mm10 (default: hg38)            
    2. **samples**:
        __The most important part of the config file is to define the samples for WES analysis.__
        Each sample is given an arbitrary name, e.g. MCF7_ER, MCF7_input, etc.  **Sample names, however, can not start with a number, and cannot contain '.', '-' (dash--use underscores instead)** (POSSIBLY others).  For each sample, define the path to the raw data file (.fastq, .fastq.gz, .bam).  For paired-end samples, simply add another line to define the path to the second pair.

    c. **setup metasheet.csv**:
    The metasheet.csv is where you group the **samples** into (defined in config.yaml) into Tunor/Normal pairs.  For WES, each of these groupings is called a **run**.

    Open metasheet.csv in Excel or in a text-editor.You will notice the first (uncommented) line is:
    `RunName,Normal,Tumor

    **RunName**- arbitraty name for the run, e.g. *MCF7_ER_run*
    **Normal**- The sample name that corresponds to Normal sample.  **It must exactly match the sample name used in config.yaml**
    **Tumor**- The tumor sample

4. editing cidc_wes/ref.yaml-
   Edit the cidc_wes/ref.yaml under the assembly for your run, i.e. hg38 or mm10.  And define the path the bwa_index file and the genome fasta file
### Running WES
1. source activate wes
2. dry run:
   A dry run helps check for any typing mistakes in the config or metasheet (or ref.yaml).  In the PROJECT directory, do this:
   snakemake -s cidc_wes/wes.snakefile -n

   If all is green, then you can move to the full run.
   If there is an error, then fix it.
3. full run:
   To do a full run, do this in the PROJECT directory:
   snakemake -s cidc_wes/wes.snakefile

   NOTE: the wes pipeline can take a long time so we recomend that you use nohup so that you can log-off and check in on the run later:
   nohup snakemake -s cidc_wes/wes.snakefile &

   NOTE: WES can use multiple cores to run.  **THIS is HIGHLY recommended**
   nohup snakemake -s cidc_wes/wes.snakefile -j [num_threads]
   where the number of threads = the number of cores you want to use.

### Appendix A: System requirements
### Appendix B: Recommended requirements
### Appendix C: Installing Chips system-wide 
###### for system administrator, those who wish to share their Chips installation
### Appendix D: Installing the mdseqpos motif finder for chips
### Appendix E: Generating static reference files for Chips
- all of the required files

- using your own files

- supporting something new

- adding to ref.yaml

  

  
