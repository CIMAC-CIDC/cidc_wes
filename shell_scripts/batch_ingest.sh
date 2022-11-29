#!/bin/bash
# Jacob Geisberg 10/2022
# A script for ingesting several samples from google bucket at once


#https://www.redhat.com/sysadmin/arguments-options-bash-scripts
Help(){
    echo "a   assay type {wes, wes_tumor_only} (required)"
    echo "d   directory containing ingestion templates (required)"
    echo "h   display help message"
}

# parse options
while getopts "ha:d:" option 
do
    case $option in
	h) Help
           exit;;
        a) ASSAY=${OPTARG};;
        d) DIRECTORY=${OPTARG};;
	\?) echo "Invalid option"
	    exit 1;;
    esac
done

# require assay and directory definitions
if [ -z "$ASSAY" ] || [ -z "$DIRECTORY" ]; then
    echo 'Missing argument'
    exit 1
fi

echo "assay: $ASSAY";
echo "directory: $DIRECTORY";

#ADD DIFFERENT COMMANDS FOR DIFFERENT ASSAYS HERE AND PASS AS VARIABLE FOR LOOP
case $ASSAY in
    wes) TXT="wes_analysis" ;;
    wes_tumor_only) TXT="wes_tumor_only_analysis";;
    *) echo "Invalid assay type"
	exit 1 ;;
esac

echo $TXT

# list of successful and failed ingestions for final output
SUCCESSES=()
FAILURES=()

#ingest each sample in specified directory. Successful templates are moved and errored templates remain 
L=$(ls $DIRECTORY)
for file in $L;
do
    echo $DIRECTORY"/"$file
    RESULT=$( yes | cidc analyses upload --analysis $TXT --xlsx $DIRECTORY"/"$file | tail -1 )

    #success
    if [[ $RESULT == *"Upload succeeded. Visit the CIDC Portal file browser to view your upload."* ]]; then
	SUCCESSES+=( $file )
	echo "success"
	if [ ! -d $DIRECTORY"_completed" ]; then
	    mkdir $DIRECTORY"_completed" 
	    echo created $DIRECTORY"_completed"
	fi
	mv $DIRECTORY"/"$file $DIRECTORY"_completed"

    #failure
    else
	FAILURES+=( $file )
	#echo "ERROR!"
	#echo $RESULT
    fi
    echo -e 
done


#final results
echo "The following templates were ingested:"
for file in "${SUCCESSES[@]}";  
do
    echo $file
done

echo -e
echo "The following templates were not ingested"
for file in "${FAILURES[@]}";
do
    echo $file
done

	    
