# Miao D, Margolis CA, Gao W, Voss MH, Li W, Martini DJ, Norton C, Bossé D,
# Wankowicz SM, Cullen D, Horak C, Wind-Rotolo M, Tracy A, Giannakis M, Hodi FS,
# Drake CG, Ball MW, Allaf ME, Snyder A, Hellmann MD, Ho T, Motzer RJ, Signoretti
# S, Kaelin WG Jr, Choueiri TK, Van Allen EM. Genomic correlates of response to
# immune checkpoint therapies in clear cell renal cell carcinoma. Science. 2018 Jan
# 4. pii: eaan5951. doi: 10.1126/science.aan5951. [Epub ahead of print] PubMed
# PMID: 29301960.
workflow captureWorkflowSNVs {
	String pairID
	String normalName
	File normalBam
	File normalBai
	String tumorName
	File tumorBam
	File tumorBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File normalPanel
	Int downsampleToCoverage
	File readgroupBlacklist
	File mutectIntervalList
	File contEstIntervalList
	File dbSNPVCF
	File dbSNPVCFIdx
	File cosmicVCF
	File picard_hapMapVCF
	File oncotatorDataSourceTarGz
	String oncotatorMode

	call ContEst_Task {
		input:
			pairID=pairID,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			tumorBam=tumorBam,
			tumorBai=tumorBai,
			normalBam=normalBam,
			normalBai=normalBai,
			contEstIntervalList=contEstIntervalList,
			picard_hapMapVCF=picard_hapMapVCF
	}
    
    call ContEstFraction_Task {
        input:
            contEst=ContEst_Task.contEst
    }

    call PrepareScatter_Task {
		input:
			targetsIntervalList=mutectIntervalList
	}

	scatter (idx in PrepareScatter_Task.scatterIndices) {
		call Mutect1_Task {
			input:
				normalName=normalName,
				normalBam=normalBam,
				normalBai=normalBai,
				tumorName=tumorName,
				tumorBam=tumorBam,
				tumorBai=tumorBai,
				refFasta=refFasta,
				refFastaIdx=refFastaIdx,
				refFastaDict=refFastaDict,
				normalPanel=normalPanel,
				downsampleToCoverage=downsampleToCoverage,
				fractionContamination=ContEstFraction_Task.contEstFraction,
				readgroupBlacklist=readgroupBlacklist,
				targetsIntervalList=PrepareScatter_Task.interval_files[idx],
				dbSNPVCF=dbSNPVCF,
				dbSNPVCFIdx=dbSNPVCFIdx,
				cosmicVCF=cosmicVCF
		}		
	}

	call GatherMutect1_Task {
		input:
			mutect1_cw=Mutect1_Task.mutect1_cw,
			mutect1_pw=Mutect1_Task.mutect1_pw,
			mutect1_cs=Mutect1_Task.mutect1_cs
	}

	call summarizeWigFile_Task {
		input:
			pairID=pairID,
			wigFile=GatherMutect1_Task.mutect1_coveragewig
	}

	call CallStatstoMAFLite_Task {
		input:
			callstats=GatherMutect1_Task.mutect1_callstats,
			pairID=pairID
	}

	call Oncotator_Task {
		input:
			MAFLITE=CallStatstoMAFLite_Task.maflite,
			oncotatorDataSourceTarGz=oncotatorDataSourceTarGz,
			oncotatorMode=oncotatorMode,
			pairID=pairID
	}

	call annotatePicardOxoQ_Task {
		input:
			pairID=pairID,
			caseBam=tumorBam,
			caseBai=tumorBai,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			dbSNPVCF=dbSNPVCF,
			dbSNPVCFIdx=dbSNPVCFIdx
	}

	call createOxoGIntervalList_Task {
		input:
			pairID=pairID,
			maf=Oncotator_Task.oncotatedMAF
	}

	call appendPicardOxoQtoMAF_Task {
		input:
			pairID=pairID,
			maf=Oncotator_Task.oncotatedMAF,
			picard_oxoQ=annotatePicardOxoQ_Task.picard_oxoQ
	}

	call createOxoGMetrics_Task {
		input:
			pairID=pairID,
			caseBam=tumorBam,
			caseBai=tumorBai,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			oxoGIntervalList=createOxoGIntervalList_Task.oxoGIntervalList
	}

	call appendOxoGtoMAF_Task {
		input:
			pairID=pairID,
			maf=appendPicardOxoQtoMAF_Task.picardOxoQMAF,
			oxoGMetrics=createOxoGMetrics_Task.oxoGMetrics
	}

	call filterOxoGArtifacts_Task {
		input:
			pairID=pairID,
			maf=appendOxoGtoMAF_Task.oxoGInfoMAF
	}

	call collectSequencingArtifacts_Task {
		input:
			pairID=pairID,
			tumorBam=tumorBam,
			tumorBai=tumorBai,
			refFasta=refFasta,
			refFastaIdx=refFastaIdx,
			refFastaDict=refFastaDict,
			dbSNPVCF=dbSNPVCF,
			dbSNPVCFIdx=dbSNPVCFIdx
	}

	call annotateOrientationBiasQ_Task {
		input:
			pairID=pairID,
			pre_adapter_detail_metrics=collectSequencingArtifacts_Task.picard_SeqArtifactMetrics
	}

	call appendFFPEQ_Task {
		input:
			pairID=pairID,
			maf=filterOxoGArtifacts_Task.oxoG3MAF,
			ffpeQ=annotateOrientationBiasQ_Task.ffpeQ
	}

	call filterFFPEArtifacts_Task {
		input:
			pairID=pairID,
			maf=appendFFPEQ_Task.ffpeQInfoMAF
	}

	call ponFilter_Task {
		input:
			pairID=pairID,
			maf=filterFFPEArtifacts_Task.ffpeBiasMAF
	}

	output {
		ContEst_Task.contEst
     	ContEst_Task.contEstTable
     	ContEst_Task.contEstBaseReport
     	ContEstFraction_Task.contEstFraction
		GatherMutect1_Task.mutect1_callstats
		GatherMutect1_Task.mutect1_powerwig
		GatherMutect1_Task.mutect1_coveragewig
		summarizeWigFile_Task.somatic_mutation_covered_bases_file_capture
        summarizeWigFile_Task.somatic_mutation_covered_bases_capture
		CallStatstoMAFLite_Task.maflite
		Oncotator_Task.oncotatedMAF
		annotatePicardOxoQ_Task.picard_oxoQ
		createOxoGIntervalList_Task.oxoGIntervalList
		appendPicardOxoQtoMAF_Task.picardOxoQMAF
		createOxoGMetrics_Task.oxoGMetrics
		appendOxoGtoMAF_Task.oxoGInfoMAF
		filterOxoGArtifacts_Task.oxoG3MAF
		filterOxoGArtifacts_Task.oxoG3PassCount
		filterOxoGArtifacts_Task.oxoG3RejectCount
		collectSequencingArtifacts_Task.picard_SeqArtifactMetrics
		annotateOrientationBiasQ_Task.ffpeQ
		appendFFPEQ_Task.ffpeQInfoMAF
		filterFFPEArtifacts_Task.ffpeBiasMAF
		filterFFPEArtifacts_Task.ffpeBiasPassCount		
		filterFFPEArtifacts_Task.ffpeBiasRejectCount
		ponFilter_Task.ponMAFAnnotated
		ponFilter_Task.ponMAFPass
		ponFilter_Task.ponBlacklist
	}
}

task ContEst_Task {
	String pairID
	File tumorBam
	File tumorBai
	File normalBam
	File normalBai
    File refFasta
    File refFastaIdx
    File refFastaDict
    File contEstIntervalList
    File picard_hapMapVCF

    command <<<
    	java -Xmx4096m	-jar /usr/GenomeAnalysisTK.jar \
    	        -T ContEst \
    	        -R ${refFasta} -I:eval ${tumorBam} -I:genotype ${normalBam} \
    	        -l INFO -pf ${picard_hapMapVCF} \
    	        -o ${pairID}.contamination.txt -br ${pairID}.contEst.baseReport.txt \
    	        -L ${contEstIntervalList} -isr INTERSECTION \
    	        --minimum_base_count 100 --trim_fraction 0.03 --beta_threshold 0.05 \
    	        --min_genotype_depth 30 --min_genotype_ratio 0.8

    	awk 'FNR == 2 {print $4}' ${pairID}.contamination.txt > ${pairID}.contEst.txt
    >>>

    runtime {
		docker: "broadinstitute/gatk3:3.8-0"
		memory: "6.5GB"
		disks: "local-disk 100 HDD"
    }

    output {
    	String contEst=read_string("${pairID}.contEst.txt") 
    	File contEstTable="${pairID}.contamination.txt"
    	File contEstBaseReport="${pairID}.contEst.baseReport.txt"
    }
}

task ContEstFraction_Task {
    String contEst

    command <<<
        python /contEstFraction.py --contEst ${contEst}
    >>>

    runtime {
        docker: "vanallenlab/contestfraction:1.0"
        memory: "1GB"
        disks: "local-disk 4 HDD"
    }

    output {
        String contEstFraction=read_string("contEstFraction.txt")
    }
}

task PrepareScatter_Task {
	File targetsIntervalList

	command <<<
		python3 /splitIntervals/splitIntervals.py --intervalHandle ${targetsIntervalList}

		# Create index for split intervals
		numIntervals="$(ls scatter.interval.* | wc -l)"
		seq 0 $((numIntervals - 1)) > indices.dat

	>>>

	output {
		Array[Int] scatterIndices=read_lines("indices.dat")
		Array[File] interval_files=glob("scatter.interval.*")
	}

	runtime {
		docker: "vanallenlab/mutect:1.1.6"
		memory: "1 GB"
	}
}

task Mutect1_Task {
	String normalName
	File normalBam
	File normalBai
	String tumorName
	File tumorBam
	File tumorBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File normalPanel
	Int downsampleToCoverage
	Float fractionContamination
	File readgroupBlacklist
	File targetsIntervalList
	File dbSNPVCF
	File dbSNPVCFIdx
	File cosmicVCF

	command <<<
		java -jar -Xmx4g /muTect-1.1.6.jar --analysis_type MuTect -L ${targetsIntervalList} \
		--normal_sample_name ${normalName} -I:normal ${normalBam} \
		--tumor_sample_name ${tumorName} -I:tumor ${tumorBam} \
		--reference_sequence ${refFasta} --normal_panel ${normalPanel} \
		--fraction_contamination ${fractionContamination} --downsample_to_coverage ${downsampleToCoverage} \
		--dbsnp ${dbSNPVCF} --cosmic ${cosmicVCF} --read_group_black_list ${readgroupBlacklist} --enable_extended_output \
		--out Mutect1.call_stats.txt --coverage_file Mutect1.coverage.wig.txt --power_file Mutect1.power.wig.txt
	>>>

	runtime {
		docker: "vanallenlab/mutect:1.1.6"
		memory: "7 GB"
		disks: "local-disk 100 HDD"
	}

	output {
		File mutect1_cs="Mutect1.call_stats.txt"
		File mutect1_pw="Mutect1.power.wig.txt"
		File mutect1_cw="Mutect1.coverage.wig.txt"
	}
}

task GatherMutect1_Task {
	Array[File] mutect1_cs
	Array[File] mutect1_pw
	Array[File] mutect1_cw

	command <<<
		MUTECT1_CW="Mutect1.coverage.wig.txt"
		cat ${sep = ' ' mutect1_cw} >> $MUTECT1_CW

		MUTECT1_PW="Mutect1.power.wig.txt"
		cat ${sep = ' ' mutect1_pw} >> $MUTECT1_PW

		MUTECT1_CS="Mutect1.call_stats.txt"
		head -2 ${mutect1_cs[0]} > $MUTECT1_CS
		cat ${sep = ' ' mutect1_cs} | grep -Pv '#' | grep -Pv '^contig' >> $MUTECT1_CS
	>>>
	
	runtime {
		docker: "broadinstitute/broadmutationcalling_beta@sha256:2285809244dabf3aefe11f4e16e91fe3d541b894f8325e4f6c9c98a8dfc81deb"
		memory: "2 GB"
	}

	output {
        	File mutect1_coveragewig="Mutect1.coverage.wig.txt"
        	File mutect1_powerwig="Mutect1.power.wig.txt"
			File mutect1_callstats="Mutect1.call_stats.txt"
	}
}

task summarizeWigFile_Task {
	String pairID
	File wigFile

	command <<<
		python /home/summarizeWigFile.py --pairId ${pairID} --wigFile ${wigFile}
	>>>

	output {
		File somatic_mutation_covered_bases_file_capture="${pairID}.somatic_coverage_summary.txt"
		String somatic_mutation_covered_bases_capture=read_string("${pairID}.somatic_coverage_summary.txt")
	}

	runtime {
		docker: "breardon/summarizewigfile:1.0"
		memory: "1 GB"
		disks: "local-disk 2 HDD"
	}
}

task CallStatstoMAFLite_Task {
	File callstats
	String pairID

	command <<<
		build=37
		mode="FSTAR"
		f_threshold=0
		triallelic_mode="REJECT"
		extra_cols="tumor_f,init_t_lod,t_lod_fstar,t_alt_count,t_ref_count,judgement"
		output="${pairID}.maf"

		perl /home/call_stats_to_maflite.pl ${callstats} $build $mode $f_threshold $triallelic_mode $output $extra_cols
	>>>

	output {
		File maflite="${pairID}.maf"
	}

	runtime {
		docker: "breardon/callstatstomaflite:1.0"
		memory: "4 GB"
		disks: "local-disk 2 HDD"
	}
}

task Oncotator_Task {
	File MAFLITE
	File oncotatorDataSourceTarGz
	String oncotatorMode
	String pairID
	File defaultConfig="gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/oncotator/tcgaMAFManualOverrides2.4.config"
	File uniProt="gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/oncotator/tx_exact_uniprot_matches.txt"

	command <<<
		tar zxvf ${oncotatorDataSourceTarGz}
		dbDir=$(basename ${oncotatorDataSourceTarGz} .tar.gz)
		mv $dbDir oncotatorDatasourceDir

		/root/oncotator_venv/bin/oncotator -i MAFLITE -o TCGAMAF \
		--db-dir=oncotatorDatasourceDir/ --tx-mode=${oncotatorMode} \
		--default_config ${defaultConfig} -v ${MAFLITE} ${pairID}.oncotated.maf hg19 \
		--log_name oncotator_firehose.log --prepend --infer-onps -c ${uniProt} --collapse-number-annotations
	>>>

	runtime {
		docker: "broadinstitute/oncotator:1.9.3.0"
		memory: "4 GB"
		disks: "local-disk 50 HDD"
	}

	output {
		File oncotatedMAF="${pairID}.oncotated.maf"
	}
}

task annotatePicardOxoQ_Task {
	String pairID
	File caseBam
	File caseBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File dbSNPVCF
	File dbSNPVCFIdx

	command <<<
		bash /1.annotatePicardOxoQ/annotatePicardOxoQ.sh -i ${pairID} -b ${caseBam} -r ${refFasta} -d ${dbSNPVCF}
	>>>
	
	runtime {
		docker: "vanallenlab/oxog_process:1.0"
		memory: "6.5GB"
		disks: "local-disk 100 HDD"
	}

	output {
		String picard_oxoQ=read_string("${pairID}.oxoQ.txt")
	}
}

task createOxoGIntervalList_Task {
	String pairID
	File maf

	command <<<
		python /2.createOxoGIntervals/createOxoGIntervals.py ${maf} ${pairID}.oxoG.interval_list
	>>>

	runtime {
		docker: "vanallenlab/oxog_process:1.0"
		memory: "4GB"
		disks: "local-disk 25 HDD"
	}

	output {
		File oxoGIntervalList="${pairID}.oxoG.interval_list"
	}
}

task appendPicardOxoQtoMAF_Task {
	String pairID
	String picard_oxoQ
	File maf

	command <<<
		bash /3.appendPicardOxoQtoMAF/AppendAnnotation2MAF.sh -i ${pairID} -m ${maf} -f "picard_oxoQ" -v ${picard_oxoQ} -o .
	>>>

	runtime {
		docker: "vanallenlab/oxog_process:1.0"
		memory: "4GB"
		disks: "local-disk 25 HDD"
	}

	output {
		File picardOxoQMAF="${pairID}.picard_oxoQ.maf.annotated"
	}
}

task createOxoGMetrics_Task {
	String pairID
	File caseBam
	File caseBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File oxoGIntervalList

	command <<<
		java -Xmx2g -jar /GenomeAnalysisTK.1.5.jar --analysis_type OxoGMetrics -R ${refFasta} -I ${caseBam} -L ${oxoGIntervalList} -o ${pairID}.oxoG.metrics.txt
	>>>

	runtime {
		docker: "vanallenlab/oxog_process:1.0"
		memory: "6.5GB"
		disks: "local-disk 100 HDD"
	}

	output {
		File oxoGMetrics="${pairID}.oxoG.metrics.txt"
	}
}

task appendOxoGtoMAF_Task {
	String pairID
	File maf
	File oxoGMetrics

	command <<<
		python /5.appendOxoGInfo/appendOxoGInfo.py --onlyAddColumnsToCopy ${oxoGMetrics} ${maf} ${pairID}.oxoGInfo.maf.annotated
	>>>

	runtime {
		docker: "vanallenlab/oxog_process:1.0"
		memory: "4GB"
		disks: "local-disk 25 HDD"
	}

	output {
		File oxoGInfoMAF="${pairID}.oxoGInfo.maf.annotated"
	}
}

task filterOxoGArtifacts_Task {
	String pairID
	File maf
	String poxoG = 0.96
	String artifactThresholdRate = 0.01
	String logThresholdRate = -1
	String oxoQ_param1 = 36
	String oxoQ_param2 = 1.5

	command <<<
		/./startFilterMAFFile /${maf} ${pairID}.oxoG3.maf.annotated ./ 0 0 ${poxoG} ${artifactThresholdRate} ${logThresholdRate} ${oxoQ_param1} ${oxoQ_param2}
	>>>

	runtime {
		docker: "vanallenlab/oxog_filter:1.0"
		memory: "4GB"
		disks: "local-disk 25 HDD"
	}

	output {
		File oxoG3MAF="${pairID}.oxoG3.maf.annotated"
		String oxoG3PassCount=read_string("${pairID}.oxoG3.maf.annotated.pass_count.txt")
		String oxoG3RejectCount=read_string("${pairID}.oxoG3.maf.annotated.reject_count.txt")
	}
}

task collectSequencingArtifacts_Task {
	String pairID
	File tumorBam
	File tumorBai
	File refFasta
	File refFastaIdx
	File refFastaDict
	File dbSNPVCF
	File dbSNPVCFIdx

	command <<<
		bash /1.collectSequencingArtifacts/CollectSequencingArtifactMetrics.sh -i ${pairID} -b ${tumorBam} -r ${refFasta} -d ${dbSNPVCF} -o .
	>>>

	runtime {
		docker: "vanallenlab/ffpe_process:1.0"
		memory: "6.5GB"
		disks: "local-disk 100 HDD"
	}

	output {
		File picard_SeqArtifactMetrics="${pairID}.pre_adapter_detail_metrics"
	}
}

task annotateOrientationBiasQ_Task {
	String pairID
	File pre_adapter_detail_metrics
	String context=".CG"
	String altBase="T"

	command <<<
		bash /2.annotateOrientationBiasQ/annotate_orientationBiasQ.sh -i ${pairID} -m ${pre_adapter_detail_metrics} -c ${context} -a ${altBase} -o .
	>>>

	runtime {
		docker: "vanallenlab/ffpe_process:1.0"
		memory: "4GB"
		disks: "local-disk 25 HDD"
	}

	output {
		String ffpeQ=read_string("${pairID}.orientation_BiasQ.txt")
	}
}

task appendFFPEQ_Task {
	String pairID
	File maf
	String annotationLabel="ffpe_Q"
	String ffpeQ

	command <<<
		bash /3.appendFFPEQ/AppendAnnotation2MAF.sh -i ${pairID} -m ${maf} -f ${annotationLabel} -v ${ffpeQ} -o .
	>>>

	runtime {
		docker: "vanallenlab/ffpe_process:1.0"
		memory: "4GB"
		disks: "local-disk 25 HDD"
	}

	output {
		File ffpeQInfoMAF="${pairID}.ffpe_Q.maf.annotated"
	}
}

task filterFFPEArtifacts_Task {
	String pairID
	File maf
	String pBias = 0.96
	String artifactThresholdRate = 0.01
	String logThresholdRate = -1
	String biasQP1 = 30
	String biasQP2 = 1.5
	String refAllele = "G"
	String artifactAllele= "A"
	String biasField = "i_ffpe"

	command <<<
		/./orientationBiasFilter /${maf} ${pairID}.ffpeBias.maf.annotated ./ 0 0 ${pBias} ${artifactThresholdRate} ${logThresholdRate} ${biasQP1} ${biasQP2} ${refAllele} ${artifactAllele} ${biasField}
	>>>

	runtime {
		docker: "vanallenlab/ffpe_filter:1.0"
		memory: "4GB"
		disks: "local-disk 25 HDD"
	}

	output {
		File ffpeBiasMAF="${pairID}.ffpeBias.maf.annotated"
		String ffpeBiasPassCount=read_string("${pairID}.ffpeBias.maf.annotated.pass_count.txt")
		String ffpeBiasRejectCount=read_string("${pairID}.ffpeBias.maf.annotated.reject_count.txt")
	}
}

task ponFilter_Task {
	String pairID
	File maf
	File pon = "gs://fc-f36b3dc8-85f7-4d7f-bc99-a4610229d66a/broadinstitute/reference/ponStewart/final_summed_tokens.hist.bin"
	String mode = "TN"
	String parameterFile = "-"
	String NMIN = 1
	String thresh = -2.5
	String WCUT = 0.5
	String CODING_ONLY = 1
	String MIN_ALT_COUNT = 3
	String Stub = ""

	command <<<
		/./maf_pon_filter /${maf} ${pon} ${pairID} ${mode} ${parameterFile} ${NMIN} ${thresh} ${WCUT} . ${CODING_ONLY} ${MIN_ALT_COUNT} ${Stub}
	>>>

	runtime {
		docker: "vanallenlab/pon_filter:1.0"
		memory: "4 GB"
		disks: "local-disk 60 HDD"
	}

	output {
		File ponMAFPass="${pairID}.pon_annotated.pass.maf"
		File ponMAFAnnotated="${pairID}.pon_annotated.maf"
		File ponBlacklist="${pairID}.pon.blacklist.txt"
	}
}
