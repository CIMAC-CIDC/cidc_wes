#WES INSTALLATION GUIDE


The WES pipeline piles together many different software tools and packages, each of which has its own set of dependancies. In fact, some of these dependancies are incompatible with each other, and in these cases, we have to create separate environments for conflicting tools. This guide aims to help the user create a working series of environments to run the entire WES pipeline.

#insert pre-reqs here (such as email about hlahd, sentieon license, mamba, conda, google cloud??, repo cloning access??, ect.) ( point to other guides for prereqs)

##Cloning WES repository and basic setup

All of the code and scripts required to run WES are in the cidc_wes repository. Our first step is to clone the repo so that we have access to it on our local system. (We like to clone our repositories in the /mnt/ssd/wes directory) To clone the repo, we can use the following command. If it works, you will see a new directory called cidc_wes.

``` shell
git clone git@bitbucket.org:plumbers/cidc_wes.git
```

#INSERT INSTRUCTIONS ABOUT HOW TO FIND A COPY OF HLAHD AND COPY IT TO AN INSTANCE!!!!

Within the cidc_wes repository, we have a setup script that significantly reduces the workload by handling all of the conda and pip installs and many other small steps needed to get the WES pipeline up and running. Unfortunately some tools still require manual steps to complete. However, we can start off the process by running the setup script. You will want to ensure that you are running the script from the directory that you would not mind copying files to. We use /mnt/ssd/wes.  If your cidc_wes repository or your hlhhd.1.4.0 is not located at /mnt/ssd/wes, you will have to open ```cidc_wes/envs/setup_wes.sh``` and change the ```HLAHD_PATH``` and ```YML_PATH``` to point to the appropriate locations. (Should probably change this) Additionally, Please note that even when using mamba, the environments may take a while to solve. Here is a command to run the setup script that will store the output in a file for you viewing in case of error.

```
nohup time bash cidc_wes/envs/setup_wes.sh > nohup_setup_wes.out &
```




## Manual steps
With our basic environments set up, we proceed to implement tools that require additional setup.

###TCELLEXTRECT
The tcellextrect module requires manually installing a few R packages that cannot be handled by conda. If not already activated, we use the following to activate the wes environment in advance of installing the R packages. This ensures that we will be able to use the packages when we use the wes environment.

```shell
source activate wes
```

Next, while the wes environment is activated, type ```R``` to activate R, and then use the following four lines of code install the four packages we need. Note that you can select any CRAN that you trust.

```R
install.packages("stringi")
install.packages("gratia")
install.packages('Rcpp')
install.packages('TcellExTRECT', repos=NULL, type='source');
```
Once you have installed all of the packages, you can quit R with ```q()```. We saved our work image when prompted. (not sure if this matters)


###XHLA
In order to implement xHLA we will have to modify some source code located inside of the bin of the XHLA conda enviroment. The first step is to navigate to the path of the bin directory of the xHLA conda environment. We can do this by listing the conda environment paths and appending ```/bin```  to the appropriate path when we change directories.

```shell
conda info --envs
# conda environments:
#
base                     /home/jacobg/miniconda3
optitype                 /home/jacobg/miniconda3/envs/optitype
sequenza                 /home/jacobg/miniconda3/envs/sequenza
vcf                      /home/jacobg/miniconda3/envs/vcf
wes                   *  /home/jacobg/miniconda3/envs/wes
xHLA                     /home/jacobg/miniconda3/envs/xHLA
```
These paths will be different on your system, and the "*" will be next to the environment you currently have activated. In our case we use the following command to change to the xHLA bin.
```shell
cd /home/jacobg/miniconda3/envs/xHLA/bin
```

Once inside the bin, We have to manually change two files to get the code to work. Please take care to make these changes in-line and do not comment anything out. The spacing in the files is poor and anything other than making direct changes may result in spacing errors. First, in ```~/miniconda3/envs/xHLA/bin/HLA/bin/align.pl```, we omit the -C 20000 argument from line 63. The new line should look like this.

```
open(IN, "diamond blastx -t . --index-mode 1 --seg no --min-score 10 --top 20 -c 1 -d $root/data/hla -q $fastq_file -f tab --quiet -o /dev/stdout |") or die $!;
```

In report.py, which is in the same directory, we have to replace ```type=file``` with ```type = argparse.FileType('r')``` in line 13 of the file. The new line should look like this:

```python
parser.add_argument('-in', dest = 'input', type = argparse.FileType('r'), required = True, help = 'Input file')
```


###VCF
Much like xHLA, the vcf environment requires a code change. First, we navigate to the bin directory of the vcf conda enviroment.

```shell
cd ~/miniconda3/envs/vcf/bin
```

From here, we need to modify lines 363-370 of a file called vcf2maf.pl, adding in a if statement. Notice both the if statement and its corresponding closing brace.

####Old code

```perl
    my $region = "$chr:" . ( $pos - 1 ) . "-" . ( $pos + length( $ref ));
    $ref_bps{$region} = $ref;
    push( @ref_regions, $region );
    $uniq_regions{$region} = 1;
    $uniq_loci{"$chr:$pos-$pos"} = 1;

}
```
####New code

```perl
    my $region = "$chr:" . ( $pos - 1 ) . "-" . ( $pos + length( $ref ));
    if ($chr ne "") {
    $ref_bps{$region} = $ref;
    push( @ref_regions, $region );
    $uniq_regions{$region} = 1;
    $uniq_loci{"$chr:$pos-$pos"} = 1;
    }
}
```
With this change complete, you have finsihed installing WES. Congratulations!!


## Notes
During the creation of the env_vars.sh files, activating and deactivating the conda environment can cause issues with the path not shutting down correctly if either file exists without the other. Please only activate the environments when neither file exists or both are completed.

If you mess up one of the steps in a way that is not easily reversible, the best approach may be to start building the offending environment from scratch. In such a case you will want to delete the environment and reconstruct it from scratch with the yaml file instead of updating the environment from yaml. For the wes enviroment the following commands would work, assuming you are in ```cidc_wes/envs``` and have the enviroment activated.
```shell
conda deactivate
mamba env remove -n wes
mamba create -n wes
mamba env update -n wes --file wes.yml
```

Snakemake is bound to snakemake>=5.30.1 because none of the channels have 5.30.1 available.

Python is bound to python==3.6.7 because later versions of python are not supported by tensorflow 2.2.2 which is a hard requirement of pvactools 2.0.7. (python versions between 3.5 and 3.8 may be ok)

Perl-bioperl, perl-dbi, are needed for both vcf2maf and perl-vcftools-vcf but these two tools are incompatible with each other due to each needing a different version of perl. This is why the vcf environment exists.

Matplotlib is bound by matplotlib<=3.0.3 since version 3.0.3 and later break the multiqc package.

Setuptools is a pip install of setuptools<58 since version 58 breaks the use_2to3 package which is needed by pvactools.

H5py 3.1.0 is installed manually because the pvactools 2.0.7 pip install will overwrite any h5py conda package and replace it with h5py 2.10.0 . However, h5py is needed for pyclone VI output. The manual install will note the incompatibility, but the install will still work.

samtools >=1.12 is used because a certain version of samtools1.8 can cause the enviroment not to solve. samtools1.12  is the version that our working environmment uses, so that is what we set as the minimum bound for samtools.
