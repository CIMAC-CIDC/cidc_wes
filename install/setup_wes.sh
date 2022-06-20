#!/bin/bash  
#Jacob Geisberg



#stores starting directory (often /mnt/ssd/wes)
export START_DIR=$PWD
#export HLAHD_PATH=$1 # add behavior for missing input?
export HLAHD_PATH=/mnt/ssd/wes/hlahd.1.4.0
export YML_PATH=/mnt/ssd/wes/cidc_wes/envs

#CLONE REPO
git clone git@bitbucket.org:plumbers/cidc_wes.git
cd cidc_wes/envs # this may not be needed if file is moved to this directory already

#BUILD ENVIRONMENT FROM YAML

mamba create -n wes
# is this needed here, use_2to3 breaks in kraken install without this but will pip work?
# or do we need to handle pip installs separately?
#pip install setuptools==57.5.0
mamba env update -n wes --file $YML_PATH/wes.yml 
# mamba env update -n wes --file wes.yml
# export CONDA_PREFIX=/michorlab/jacobg/miniconda3
# export CONDA_ROOT=/michorlab/jacobg/miniconda3
# export PATH=/michorlab/jacobg/miniconda3/bin:$PATH
source activate wes
pip install setuptools==57.5.0
pip install pvactools==2.0.7 vcf-annotation-tools 
conda deactivate


mamba create -n xHLA
mamba env update -n xHLA --file $YML_PATH/xHLA.yml

mamba create -n optitype
mamba env update -n optitype --file $YML_PATH/optitype.yml

mamba create -n sequenza
mamba env update -n sequenza --file $YML_PATH/sequenza.yml

mamba create -n vcf
mamba env update -n vcf --file $YML_PATH/vcf.yml


#MSISENSOR
source activate wes # added to get CONDA_PREFIX variable
cd $START_DIR
git clone https://github.com/niu-lab/msisensor2.git
cd msisensor2
cp msisensor2 $CONDA_PREFIX/bin


#HLAHD
# maybe make path to hlahd an argument? Otherwise we can't ensure successful copy
cd $CONDA_PREFIX/bin
#cp -r /home/taing/miniconda3/envs/wes/bin/hlahd.1.4.0 . # replace path with the path to hlahd install 
cp -r $HLAHD_PATH .
ln -s hlahd.1.4.0 hlahd
#setting path variable done at the end of the wes section


#TCELLEXTRECT
cd $START_DIR
git clone https://github.com/McGranahanLab/TcellExTRECT.git
#manual steps to activate wes and install stringi, gratia, Rcpp and tcellextrect R packages 


#CLONALITY
pip install git+https://github.com/Roth-Lab/pyclone-vi.git
pip install h5py==3.1.0

#NEOANTIGEN
sudo apt-get install gcc-multilib g++-multilib # should this be done before conda env?

# ADD WES PATHING HERE
export WES_PATH=$CONDA_PREFIX
conda deactivate # wes must be deactivated for this part!                                                                 
cd $WES_PATH/etc/conda/activate.d
# echo """export PREPATH=\$PATH                                                                                              
# export PATH=\$PATH:~/miniconda3/envs/wes/bin/hlahd                                                                         
# """ > env_vars.sh
echo """export PREPATH=\$PATH
export PATH=\$PATH:$WES_PATH/bin/hlahd
""" > env_vars.sh
cd $WES_PATH/etc/conda/deactivate.d
echo """export PATH=\$PREPATH                                                                                              
unset PREPATH                                                                                                             
""" > env_vars.sh

#XHLA (section should be tested)
cd $START_DIR
git clone https://github.com/humanlongevity/HLA.git
source activate xHLA #needed for CONDA_PREFIX and for diamond
cd $CONDA_PREFIX/bin
cp -R $START_DIR/HLA .
cd HLA/data
diamond makedb --in hla.faa -d hla
#manual file changes to aligner.pl and report.py in the ~/miniconda3/envs/xHLA/bin/HLA/bin directory


export XHLA_PATH=$CONDA_PREFIX
conda deactivate
cd $XHLA_PATH/etc/conda/activate.d
# echo """export PREPATH=\$PATH                                                                                              
# export PATH=\$PATH:~/miniconda3/envs/xHLA/bin/HLA/bin                                                                      
# """ > env_vars.sh
echo """export PREPATH=\$PATH
export PATH=\$PATH:$XHLA_PATH/bin/HLA/bin 
""" > env_vars.sh
cd $XHLA_PATH/etc/conda/deactivate.d
echo """export PATH=\$PREPATH 
unset PREPATH                                                                                    
""" > env_vars.sh

#SEQUENZA
conda activate sequenza
cd $CONDA_PREFIX/bin 
#ln -s ~/miniconda3/envs/wes/bin/x86_64-conda_cos6-linux-gnu-cc
ln -s $WES_PATH/bin/x86_64-conda_cos6-linux-gnu-cc #untested on google cloud but should replace above
conda deactivate

#VCF
#manual file changes to ~/miniconda3/envs/vcf/bin/vcf2maf.pl

