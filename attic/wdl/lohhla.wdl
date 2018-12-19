workflow run_lohhla {
	call lohhla
}

task lohhla {
	String patientID
	File allBAM
	File hla
	File HLAfastaLoc
	File CopyNumLoc
	Int? minCoverageFilter = 10
	Boolean? mappingStep = true
	Boolean? fishingStep = true
	Boolean? cleanUp = false
	Boolean? plottingStep = true
	File? normalBam
	
	# system
	Int? memory = 16
	Int? disk_space = 50
	Int? num_cpu = 1
	Int? num_preempt = 0
	Int? boot_disk_gb = 30
	String? docker = "gcr.io/cidc-dfci/lohhla:2.0"
	command <<<
		#set -euo pipefail # bash strict mode
		mkdir BAM
		tar -xvvf ${allBAM} -C BAM --strip-components=1
		mkdir lohhlaOut
		if [ ! -f "${normalBam}"  ];then
			normalFile=FALSE
		else
			mv ${normalBam} BAM/
			normalFile=$PWD/BAM/$(basename "${normalBam}")
		fi

		Rscript /opt/lohhla/LOHHLAscript.R --patientId ${patientID} --outputDir $PWD/lohhlaOut/ --BAMDir $PWD/BAM/ --hlaPath ${hla} --normalBAMfile $normalFile \
--HLAfastaLoc ${HLAfastaLoc} \
--CopyNumLoc ${CopyNumLoc} \
--mappingStep ${mappingStep} \
--minCoverageFilter ${minCoverageFilter} \
--fishingStep ${fishingStep} \
--cleanUp ${cleanUp} \
--gatkDir /opt/picard-tools-1.119/ \
--novoDir /opt/novocraft/ \
--HLAexonLoc /opt/lohhla/data/hla.dat \
--plottingStep ${plottingStep} 
	>>>
	
	output {
	File HLAlossPrediction="lohhlaOut/${patientID}.${minCoverageFilter}.DNA.HLAlossPrediction_CI.xls"
	}
	
	meta {
	maintainer:"Jingxin Fu"
	citation: "McGranahan et al., Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution, Cell (2017)"
	}
	
	runtime {
	docker:"${docker}"
	zones: "us-east1-b us-east1-c us-east1-d"
	memory: "${memory} GB"
        #slurm_memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${num_preempt}"
        bootDiskSizeGb: "${boot_disk_gb}"
	}
}
