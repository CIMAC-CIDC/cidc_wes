#!/usr/bin/env python

from genson import SchemaBuilder

###############################################################################
# Constants
###############################################################################
_int = 777
_float = 7.0
_string = "TGBTG"
_int_arr = [_int,_int,_int]

###############################################################################
# SAMPLE level definition
###############################################################################
alignment = {'total_reads': _int,
             'mapped_reads': _int,
             'dedup_reads': _int,
             'gc_content': _int_arr, #NOTE: could be floats
             'quality_score': _int_arr} #NOTE: could be floats

coverage = {'total_reads': _int,
            'mean_depth': _float,
            "q1_depth": _float,
            "median_depth": _float,
            "q3_depth": _float,
            "percent_bases_gt_50": _float}

hla = { "A-1": _string,
        "A-2": _string,
        "B-1": _string,
        "B-2": _string,
        "C-1": _string,
        "C-2": _string,
        "DPB1-1": _string, #OPTIONAL
        "DPB1-2": _string, #OPTIONAL
        "DQB1-1": _string, #OPTIONAL
        "DQB1-2": _string, #OPTIONAL
        "DRB1-1": _string, #OPTIONAL
        "DRB1-2": _string} #OPTIONAL

sample = {'id': _string,
          'alignment': alignment,
          'coverage': coverage,
          'hla': hla}

###############################################################################
# END SAMPLE level definition
###############################################################################

###############################################################################
# SOMATIC level definition
###############################################################################

mutation_results = {"total": _int,
                    "snp": _int,
                    "insertion": _int,
                    "deletion": _int}

transition_row = {"A": _int,
                  "C": _int,
                  "G": _int,
                  "T": _int}


somatic_results = {'mutation_summary': mutation_results,
                   'functional_summary': mutation_results,
                   'trinucleotide_matrix': _int_arr,
                   'transition_matrix': {'A': transition_row,
                                         'C': transition_row,
                                         'G': transition_row,
                                         'T': transition_row}}
###############################################################################
# END SOMATIC level definition
###############################################################################

###############################################################################
# RUN level definition
###############################################################################

neoantigen_row = {"Gene": _string,
                  "EnsemblID": _string,
                  "HLA": _string,
                  "Peptide_Sequence": _string,
                  "Read_Depth": _float,
                  "DNA_VAF": _float,
                  "Method": _string,
                  "Score": _float,
                  "WT_Score": _float,
                  "Fold_Change": _float}

run = {'id': _string,
       'tumor': sample,
       'normal': sample,
       'somatic': somatic_results,
       'neoantigen': [neoantigen_row, neoantigen_row, neoantigen_row],
       "cnv": _string, #TODO
       "clonality": _float, #Q: can clonality be represented by a single #?
       "germline": _float,
       "purity": _float,
       }

###############################################################################
# END RUN level definition
###############################################################################

builder = SchemaBuilder(schema_uri="http://json-schema.org/draft-07/schema#")
builder.add_object(run)
print(builder.to_json(indent=3))
