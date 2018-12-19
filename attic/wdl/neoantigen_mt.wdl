# ******************************************************** #
#   WDL Script for Neoantigen Pipeline with NetMHCpan-4.0  #
#   With the Multiple Threads Improvement                  #
# ******************************************************** #

workflow Neoantigen {	

	# RUNTIME INPUT PARAMS
	# non-negative interger value for preemptible 0 means not preemptible, 
	# otherwise 1,2,... is the max number of pre-emptible tries
	Int preemptible
	
	# WORKFLOW INPUT PARAMS
	# directory for hg19 - /firecloud-tcga-open-access/tutorial/reference/annotation.db.ucsc.hg19.tar
	File hg19DBTarBall
	# SNVs annotated MAF file
	File mutFile 
	# INDELs annotated MAF file
	File indelFile
	# HLA file 
	File hlaFile
	# identifier of the individual whose mutFile and indelFile are processed by this script
	String id 
	# ranking scheme to use, default = 2
	Int? rankVersion
	# all or std
	String binderThreshold
	# keep or remove HLA-C predictions (0==remove, 1=keep)
	Int hlaCStatus

	call get_peptides {
		input: 
			hg19DBTarBall = hg19DBTarBall,
			mutFile = mutFile,
			indelFile = indelFile,
			hlaFile = hlaFile,
			id = id,			
			preemptible=preemptible
	}

    call split_files {
        input:
		    outFilePep17 = get_peptides.outFilePep17,
		    outFilePep19 = get_peptides.outFilePep19,
		    outFileIndelPep17 = get_peptides.outFileIndelPep17,
		    outFileIndelPep19 = get_peptides.outFileIndelPep19
    }

    scatter (pepFile in split_files.pepFiles) {
        call predict_peptides {
            input: 
                hlaForNetMHC = get_peptides.hlaForNetMHC,
                pepFile = pepFile,
                name = basename(pepFile),
                binderThreshold = binderThreshold,
                preemptible=preemptible
        }
    }

    call merge_results {
        input:
            binders = predict_peptides.binders,
            outFullPeptide = get_peptides.outFullPeptide,
            outIndelFullPeptide = get_peptides.outIndelFullPeptide,
            outFilePep17 = get_peptides.outFilePep17,
            id = id,
	        rankVersion = rankVersion,
            binderThreshold = binderThreshold,
            hlaCStatus = hlaCStatus
    }

	output {			
        File hlaForNetMHC = get_peptides.hlaForNetMHC
		File outFullPeptide = get_peptides.outFullPeptide
		File outFilePep17 = get_peptides.outFilePep17
		File outFilePep19 = get_peptides.outFilePep19
		File outFilePep33 = get_peptides.outFilePep33
		File outFilePep81 = get_peptides.outFilePep81
		File outIndelFullPeptide = get_peptides.outIndelFullPeptide
		File outFileIndelPep17 = get_peptides.outFileIndelPep17
		File outFileIndelPep19 = get_peptides.outFileIndelPep19
		File outFileIndelPep33 = get_peptides.outFileIndelPep33
		File outFileIndelPep81 = get_peptides.outFileIndelPep81

		File combinedFullpeptideFile = merge_results.combinedFullpeptideFile
		File combinedBindersAnnotFile = merge_results.combinedBindersAnnotFile
		File neoorfFile = merge_results.neoorfFile
		File combinedAllBindersCleanFile = merge_results.combinedAllBindersCleanFile
		File outFilterInconsistentFile = merge_results.outFilterInconsistentFile
	}

}

task get_peptides {

	# RUNTIME INPUT PARAMS
	Int? preemptible = 0
    String? docker = "gcr.io/cidc-biofx/neoantigen:v4"

	# TASK INPUT PARAMS
	File hg19DBTarBall
	File mutFile
	File indelFile
	File hlaFile
	String id	

	# COMPUTE DISK SIZE
	Int diskGB = ceil(size(hg19DBTarBall, 'G') * 4 + size(mutFile, 'G') + size(indelFile, 'G') + 20) 

	command <<<

        # format HLA name
        cat ${hlaFile} | awk '{printf("%s\n%s\n",$2,$3)}' | uniq | sed 's/hla_a/HLA-A/g' | sed 's/hla_b/HLA-B/g' | sed 's/hla_c/HLA-C/g' | sed 's/\(HLA-[ABC]\)_\(..\)_\(..\).*/\1\2:\3/g' > ${id}.hla.txt

        # clean up mutation file and limit computed cases
        cat ${mutFile} | egrep '^#version|^## Oncotator|^Hugo_Symbol|Missense_Mutation\sSNP' > mutFile.maf
        cat ${indelFile} | egrep '^#version|^## Oncotator|^Hugo_Symbol|Frame_Shift_Del\sDEL|Frame_Shift_Ins\sINS' > indelFile.maf

		set -euxo pipefail		

        HG19_DB_DIR_NAME=`pwd`/tempDir
        if [[ -d /mnt/ramdisk ]]; then
            HG19_DB_DIR_NAME=/mnt/ramdisk/tempDir
        fi

        mkdir $HG19_DB_DIR_NAME
        tar zxf ${hg19DBTarBall} -C $HG19_DB_DIR_NAME/ --strip-components 7

		# modified get_mut_protein_gencode_v1.pl		
		/root/neoantigen/src/get_mut_protein_gencode_v1.pl mutFile.maf infile transcript_gencode_v19 ${id} $HG19_DB_DIR_NAME > ${id}.mutFile.out
		/root/neoantigen/src/get_mut_protein_gencode_v1.pl indelFile.maf infile transcript_gencode_v19 ${id}.indel $HG19_DB_DIR_NAME > ${id}.indelFile.out

	>>>

	runtime {
		docker: "${docker}"
		slurm_docker: "${docker}"
		memory: "16 GB"
		slurm_memory: "16G"
		disks: "local-disk ${diskGB} HDD"
		preemptible: "${preemptible}"
        zones: "us-east1-b us-east1-c us-east1-d"
        mount_tmpfs: "/mnt/ramdisk"
	}

	output {			
        File hlaForNetMHC = "${id}.hla.txt"
		File outFullPeptide = "${id}.fullpeptide.fa" 
		File outFilePep17 = "${id}.pep17.fa" 
		File outFilePep19 = "${id}.pep19.fa" 
		File outFilePep33 = "${id}.pep33.fa" 
		File outFilePep81 = "${id}.pep81.fa"
		File outIndelFullPeptide = "${id}.indel.fullpeptide.fa" 
		File outFileIndelPep17 = "${id}.indel.pep17.fa" 
		File outFileIndelPep19 = "${id}.indel.pep19.fa" 
		File outFileIndelPep33 = "${id}.indel.pep33.fa" 
		File outFileIndelPep81 = "${id}.indel.pep81.fa" 		
	}
}


task split_files {

	File outFilePep17
	File outFilePep19
	File outFileIndelPep17
	File outFileIndelPep19

    Int? length=800 ## must be the multiply of 8
    String? docker = "gcr.io/cidc-biofx/neoantigen:v4"

    command <<<

        split -l ${length} ${outFilePep17} pep17
        split -l ${length} ${outFilePep19} pep19
        split -l ${length} ${outFileIndelPep17} pepin17
        split -l ${length} ${outFileIndelPep19} pepin19

        for file in pep*
        do
            mv "$file" "$file.fa"
        done
        
    >>>

	runtime {
		docker: "${docker}"
		slurm_docker: "${docker}"
        zones: "us-east1-b us-east1-c us-east1-d"
		disks: "local-disk 60 HDD"
		
	}

    output {
        Array[File] pepFiles = glob("*.fa")
    }
}


task predict_peptides {

	# RUNTIME INPUT PARAMS
	Int? preemptible = 0
    String? docker = "gcr.io/cidc-biofx/neoantigen:v4"

	# TASK INPUT PARAMS
    File hlaForNetMHC
	File pepFile
	String name
	String binderThreshold		

	command <<<
		
        pepFile="${name}"
        cp ${pepFile} $pepFile

		/root/neoantigen/src/get_netmhc_results_cmds_netmhc4.pl $pepFile 9 ${hlaForNetMHC} /root/neoantigen/all_netmhc_variants.txt mut 
		sh $pepFile".nmhc.cmd" # Calling NetMHCpan-4.0 - eluted ligand likelihood
		sh $pepFile".BA.nmhc.cmd" # Calling NetMHCpan-4.0 - binding affinity
		/root/neoantigen/src/get_netmhc_results_post_run_netmhc4.pl $pepFile 9 ${hlaForNetMHC} /root/neoantigen/all_netmhc_variants.txt mut ${binderThreshold}

		/root/neoantigen/src/get_netmhc_results_cmds_netmhc4.pl $pepFile 10 ${hlaForNetMHC} /root/neoantigen/all_netmhc_variants.txt mut 
		sh $pepFile".nmhc.cmd" # Calling NetMHCpan-4.0 - eluted ligand likelihood
		sh $pepFile".BA.nmhc.cmd" # Calling NetMHCpan-4.0 - binding affinity
		/root/neoantigen/src/get_netmhc_results_post_run_netmhc4.pl $pepFile 10 ${hlaForNetMHC} /root/neoantigen/all_netmhc_variants.txt mut ${binderThreshold}

	>>>

	runtime {
		docker: "${docker}"
		slurm_docker: "${docker}"
		memory: "16 GB"
		slurm_memory: "16G"
        cpu: "2"
		preemptible: "${preemptible}"
        zones: "us-east1-b us-east1-c us-east1-d"
		disks: "local-disk 60 HDD"
	}

	output {			
        File binders = "${name}.${binderThreshold}.binders.txt"
	}
}
	

task merge_results {

	# RUNTIME INPUT PARAMS
	Int? preemptible = 0
    String? docker = "gcr.io/cidc-biofx/neoantigen:v4"

	# TASK INPUT PARAMS
	String id
    Array[File] binders
    File outFullPeptide
    File outIndelFullPeptide
    File outFilePep17

	Int? rankVersion
	String binderThreshold		
	Int hlaCStatus

	command <<<

		echo "#### aggregate 9 and 10 mer binders ####"

        for binder in ${sep=' ' binders}
        do
            if [ ! -f ${id}.combined.${binderThreshold}.binders.txt ]
            then 
                head -1 $binder > ${id}.combined.${binderThreshold}.binders.txt
            fi
            tail -n+2 $binder  >> ${id}.combined.${binderThreshold}.binders.txt	 
        done

		cat ${outFullPeptide} ${outIndelFullPeptide} > ${id}.combined.fullpeptide.fa

		echo "#### process_binders ####"
		/root/neoantigen/src/process_binders_netmhc4.pl ${id}.combined.${binderThreshold}.binders.txt ${id}.combined.fullpeptide.fa  /root/neoantigen/wt.peptide.9mers.txt /root/neoantigen/wt.peptide.10mers.txt mut ${default=2 rankVersion} ${hlaCStatus} > ${id}.${binderThreshold}.process_binders.out

		echo "#### get_neoorf_table ####"
		/root/neoantigen/src/get_neoorf_table.pl ${id}.combined.fullpeptide.fa /root/neoantigen/wt.peptide.9mers.txt /root/neoantigen/wt.peptide.10mers.txt
		
		echo "#### filter inconsistent ####"
		/root/neoantigen/src/filter_inconsistent.pl ${outFilePep17} ${id}.combined.${binderThreshold}.binders.txt.annot.txt > ${id}.inconsistent.mutations.txt
	
	>>>

	runtime {
		docker: "${docker}"
		slurm_docker: "${docker}"
		memory: "32 GB"
		slurm_memory: "32G"
		preemptible: "${preemptible}"
        zones: "us-east1-b us-east1-c us-east1-d"
		disks: "local-disk 60 HDD"
	}

	output {
		File combinedFullpeptideFile = "${id}.combined.fullpeptide.fa" 			
		File combinedBindersAnnotFile = "${id}.combined.${binderThreshold}.binders.txt.annot.txt" 
		File neoorfFile = "${id}.combined.neoorf.table.txt"
		File combinedAllBindersCleanFile = "${id}.combined.${binderThreshold}.binders.txt.annot.txt.clean.txt"
		File outFilterInconsistentFile = "${id}.inconsistent.mutations.txt"	
	}
}
