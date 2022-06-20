#WES INSTALLATION GUIDE


The WES pipeline piles together many different software tools and packages, each of which has its own set of dependancies. In fact, some of these dependancies are incompatible with each other, and in these cases, we have to create separate environments for conflicting tools. This guide aims to help the user create a working series of environments to run the entire WES pipeline.

#insert pre-reqs here (such as email about hlahd, sentieon license, mamba, conda, google cloud??, repo cloning access??, ect.) ( point to other guides for prereqs)

##Cloning WES repository and basic setup

All of the code and scripts required to run WES are in the cidc_wes repository. Our first step is to clone the repo so that we have access to it on our local system. (We like to clone our repositories in the /mnt/ssd/wes directory) To clone the repo, we can use the following command. If it works, you will see a new directory called cidc_wes.

``` shell
git clone git@bitbucket.org:plumbers/cidc_wes.git
```

Within the cidc_wes repository, there are several yaml files that can be used to create the conda environments that the WES pipeline needs to run. Unfortunately, not all tools are available via conda, and some tools require additional steps to get to a working state. However, we can start off the process by navigating to the cidc_wes/envs directory and creating the main wes environment using the following commands. Please note that even when using mamba, the environment may take a while to solve.

```shell
cd cidc_wes/envs
mamba create -n wes
mamba env update -n wes --file wes.yml
```

The first command puts us in the cidc_wes/envs folder where the yml files are located. The second and third commands create and update the environment. While it is possible to combine the latter two commands with ```mamba env create -f wes.yml```, we do not recommend this as it still reverts to conda for the installation of packages.

Next, the wes environment has some packages that are installed by pip. To install them, we activate the environment, run the installation commands, and then deactivate the environment. Critically, setuptools is installed first as it is needed fot the latter install.

```shell
source activate wes
pip install setuptools==57.5.0
pip install pvactools==2.0.7 vcf-annotation-tools 
conda deactivate
```

Once the basic setup for wes environment is finished, we create the baseline for our other conda environments, xHLA, optitype, sequenza and vcf. While there may be other yml files present in the envs folder, these are the only conda environments needed to run the WES pipeline.

```shell
mamba create -n xHLA
mamba env update -n xHLA --file xHLA.yml

mamba create -n optitype
mamba env update -n optitype --file optitype.yml

mamba create -n sequenza
mamba env update -n sequenza --file sequenza.yml

mamba create -n vcf
mamba env update -n vcf --file vcf.yml
```



## WES ENVIRONMENT
With our basic environments set up, we proceed to implement tools that require additional setup, starting with those that will be run from the main wes environment. The steps are generally broken down by which module they pertain to, but some modules require multiple conda environments, and some conda environments are needed for multiple modules.  While a few of these steps *could* be skipped depending on the user's individual needs, we strongly recommend completing all of them **in the order they appear here**. This is the only way to be assured that the WES pipeline will run properly.



<!--
###SOMATIC
First, we install mamba install perl-vcftools, which needs to be done while the wes environment is active. If it is not already active, any conda environment can be activated with ``` source activate <name> ``` where <name> is the sname of the environment you would like to activate. Conversely any environment can be deactivated with ```conda deactivate <name>```.

<!-- In our case we activate the wes environment then install perl-vcftools with mamba.

```shell
source activate wes
mamba install -c bioconda perl-vcftools-vcf -n wes
```

Next, we have to slightly modify the source code of vcf2maf in order to make it work. The first step is to navigate to the path of the bin directory of the wes conda environment. We can do this by listing the conda environment paths and appending ```/bin```  to the appropriate path when we change directories.

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

cd /Users/jacobg/miniconda3/envs/wes/bin
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
-->
###MSISENSOR
For msisensor, we start by cloning the msisensor2 repository from github, and navigating into it.

```shell
git clone https://github.com/niu-lab/msisensor2.git
cd msisensor2
```

Once inside the new msisensor2 directory, we have to copy a script named msisensor2 into the bin directory of the wes conda environment. We can identify the path to the wes bin by listing the conda environment paths and appending ```/bin```  to the appropriate path.

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
In our case we have ```/home/jacobg/miniconda3/envs/wes/bin``` which we can replace with ```~/miniconda3/envs/wes/bin``` . Now that we have the path identified, we can copy the msisensor2 script over to it.

```shell
cp msisensor2 ~/miniconda3/envs/wes/bin
```


###HLAHD (note that proper install/email may be covered elsewhere in doc)


This section assumes you have been able to acquire and install a copy of HLAHD on your system. We copy the hlahd.1.4.0 into the wes bin. Then, we use the ln command to create a symbolic link to hlahd.1.4.0 called hlahd.

```shell
cd ~/miniconda3/envs/wes/bin
cp -r /home/taing/miniconda3/envs/wes/bin/hlahd.1.4.0 .
ln -s hlahd.1.4.0 hlahd
```

We also have to create files to add and remove the symbolic link from our path variable. This will be covered in our finishing touches section.

<!-- Next, we have to add the symbolic link that we created to our PATH. We would like to do this automatically every time wes conda enviroment is activated so we need to navigate to the etc/conda/activate.d folder of our wes enviroment. Once there, we create a file called env_vars.sh. (not sure the ~ will work here...)

```shell
cd ~/miniconda3/envs/wes/etc/conda/activate.d
emacs env_vars.sh
```

Inside the file, add the following two lines of bash code.

``` bash
export PREPATH=$PATH
export PATH=$PATH:~/miniconda3/envs/wes/bin/hlahd
```

Finally, we need to remove hlahd from our path if we deactivate the environment. To do this, we change into the conda/deactivate.d Folder and add the following two lines of code to a new file called env_vars.sh

```bash

``` -->


###TCELLEXTRECT
To use TcellExTRECT, we have to clone the TcellExTRECT repository and manually install some R packages. To install the repository we use:

```shell
cd /mnt/ssd/wes
git clone https://github.com/McGranahanLab/TcellExTRECT.git
```

If not already activated, we use the following to activate the wes environment in advance of installing the R packages.

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
You can quit R with ```q()```. We saved our work image when prompted. (not sure if this matters)


### CLONALITY
The clonality module runs two progams, Pyclone VI and Seqeunza. Sequenza runs in its own conda environment, so we will cover that install later. For Pyclone VI, we need to install its repository while wes is active. we also need to manually install h5py==3.1.0 . Note that the h5py install is done manually because of pvactools requiring tensorflow 2.2.2 which requires h5py 2.10 . As a result, you may see an error here, but install should still work. (replace with reference to things that can go wrong)

```shell
pip install git+https://github.com/Roth-Lab/pyclone-vi.git
pip install h5py==3.1.0
```

You can test that the h5py version is correct with ```conda list | grep "h5py"``` .


###NEOANTIGEN
All that is needed for Neoantigen is the following command. (should this be done before conda env??)

```shell
sudo apt-get install gcc-multilib g++-multilib
```



##xHLA
Once we have already created the xHLA enviroment, the first step is to clone the xHLA repository.

```shell
git clone https://github.com/humanlongevity/HLA.git
```

Next we copy the resulting HLA directory to the bin directory of the xHLA environment.

```shell
cd ~/miniconda3/envs/xHLA/bin
cp -r /mnt/ssd/wes/HLA .
```

From inside xHLA/bin/HLA/data , we create the diamond db. Note that we have to have the xHLA environment activated for this step. (I always deactivate one environment before activating another, is this necessary)

```shell
source activate xHLA
diamond makedb --in hla.faa -d hla
```

<!-- Much as we did with hlahd, we have to add xHLA to our path when we activate the xHLA enviroment. To do this we go to the etc/conda/activate.d directory and add the folling to a new env_vars.sh file (replace ~ ??)
```bash
export PREPATH=$PATH
export PATH=$PATH:~/miniconda3/envs/xHLA/bin/HLA/bin
``` -->

Finally, We have to manually change two files to get the code to work. Please take care to make these changes in-line and do not comment anything out. The spacing in the files is poor and anything other than making direct changes may result in spacing errors. First, in ```~/miniconda3/envs/xHLA/bin/HLA/bin/align.pl```, we omit the -C 20000 argument from line 63. The new line should look like this.

```
open(IN, "diamond blastx -t . --index-mode 1 --seg no --min-score 10 --top 20 -c 1 -d $root/data/hla -q $fastq_file -f tab --quiet -o /dev/stdout |") or die $!;
```

In report.py, which is in the same directory, we have to replace ```type=file``` with ```type = argparse.FileType('r')``` in line 13 of the file. The new line should look like this:

```python
parser.add_argument('-in', dest = 'input', type = argparse.FileType('r'), required = True, help = 'Input file')
```


##SEQUENZA
The sequenza environment is pretty straightforward to set up, as it only requires us to create a link between the the gcc install in our wes enviroment and the sequenza bin.

```shell
cd ~/miniconda3/envs/sequenza/bin
ln -s ~/miniconda3/envs/wes/bin/x86_64-conda_cos6-linux-gnu-cc
```

##VCF
First, we navigate to the bin directory of the vcf conda enviroment
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
#INSERT OPTITYPE CHANGES HERE!!!!!!!!!!!!!!!!!!!


#Finishing Touches
Our wes and xhla environments require that we add variables to our path when we activate them, and remove those variables upon deactivation. In both cases we create files at ~/miniconda3/envs/<env_name>/etc/conda/activate.d/env_vars.sh and ~/miniconda3/envs/<env_name>/etc/conda/deactivate.d/env_vars.sh . These handle activation and deactivation respectively. To avoid annoying errors later, we will want to deactivate conda and keep it deactivated until both files are in place. The command for deactivating conda is as follows.

```shell
conda deactivate
```


####wes

For wes, we change to the appropriate directory and create env_vars.sh. (do ~ work here?)

```shell
cd ~/miniconda3/envs/wes/etc/conda/activate.d
emacs env_vars.sh
```

Inside the file, add the following two lines of bash code. These are needed to ensure that HLAHD works correctly.

``` bash
export PREPATH=$PATH
export PATH=$PATH:~/miniconda3/envs/wes/bin/hlahd
```

Next, we change to the deactivate.d directory and create another file called env_vars.sh

```shell
cd ~/miniconda3/envs/wes/etc/conda/deactivate.d
emacs env_vars.sh
```

We add the following two lines of code to deactivate the environment

```bash
export PATH=$PREPATH
unset PREPATH
```
####xHLA
For xhla we follow a similar process to above.

```shell
cd ~/miniconda3/envs/xHLA/etc/conda/activate.d
emacs env_vars.sh
```
Inside the file add:

```bash
export PREPATH=$PATH
export PATH=$PATH:~/miniconda3/envs/xHLA/bin/HLA/bin
```
Our deactivation file is the same as the one used for the wes enviroment

```shell
cd ~/miniconda3/envs/xHLA/etc/conda/deactivate.d
emacs env_vars.sh
```

```bash
export PATH=$PREPATH
unset PREPATH
```


<!-- deactivate here
deactivate and reactivate!!! -->

<!-- cd ~/miniconda3/envs/<env of choice>/etc/conda/deactivate.d
emacs env_vars.sh

```bash
export PATH=$PREPATH
unset PREPATH
``` -->

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

In optitype.yml, pysam and and razers3 are set to pysam==0.15.2 and razers3==3.5.3 becuase pysam 0.11.2.2 and razors3 3.5.8 give "Warning: PySam not available on the system. Falling back to primitive SAM parsing."	which slows the	optitype module down by multiple hours in some cases.