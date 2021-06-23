/* Len Taing 2021 (TGBTG) 
   Script to decode the somatic.filtered_maf_file from base64 encoding 
   and also provide maf information for other js components
*/

/* mafStringToDf- 
   takes a string representing the contents of an entire maf file and 
   returns a danfo DataFrame containing the maf file information

   This fn will also select out the key fields that are used in the cohort
   report.

   ref: https://stackoverflow.com/questions/16177037/how-to-extract-information-in-a-tsv-file-and-save-it-in-an-array-in-javascript/16177134
*/
/* OBSOLETE and inefficient!
function mafStringToDf(mafString) {
    let lines = mafString.split("\n"); //split maf files into llines
    let mat = []; //matrix of maf contents
    let hdr = null;
    for (var i=0; i<lines.length; i++) {
        if (lines[i].substring(0,1) != "#") { //skip comments
	    //Check for special case, first read in is hdr
	    if (hdr == null) {
		hdr = lines[i].split('\t');
	    } else {
		mat.push(lines[i].split('\t'));
	    }
	}
    }

    //MAKE df object
    let df = new dfd.DataFrame(mat, {'columns': hdr});
    //Select only relevant cols here
    const _fields = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 't_ref_count', 't_alt_count', 'HGVSc', 'HGVSp', 'Transcript_ID'];
    subset_df = df.loc({columns: _fields})

    return subset_df;
}
*/

//Convert all wes_data[i]['somatic']['filtered_maf_file'] into a danfo Df
//Assumes that each wes_data[i]['somatic']['filtered_maf_file'] object is
// {'hdr' : [ ... ], 'mat': [[...], [ ...], ]}
for (var i = 0; i < wes_data.length; i++) {
    let tmp = wes_data[i]['somatic']['filtered_maf_file'];
    let df = new dfd.DataFrame(tmp['mat'], {'columns': tmp['hdr']});
    wes_data[i]['somatic']['maf'] = df;
    //delete filtered_maf_file
    delete wes_data[i]['somatic']['filtered_maf_file'];
}

//EXAMPLES:
//Print the Chromosome cols from the 1st wes run
//wes_data[0].somatic.maf["Chromosome"].print();
