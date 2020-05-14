#!/usr/bin/env python3
"""Script to generate analysis/report/wes_level1.html"""

import os
import sys

import re
import subprocess
import yaml
import json

from optparse import OptionParser
import jinja2
import pandas as pd

from report_level1 import getFileName
from report_level1 import getRuns
from clonality_calcClonality import calcClonality

################ THIS fn is copied from neoantigen.snakefile ##################
def parseHLA(config, hla_files):
    """Given an optitypes '_results.tsv' file; parses the HLA A, B, C
    and returns these as a comma-separated string (for pvacseq) input

    NOTE: cureently the optitype results.tsv looks somthing like this:
    	A1	A2	B1	B2	C1	C2	Reads	Objective
    0					C*06:04	C*06:04	4.0	3.99
    **So were' going to parse cols 1-6 and return that"""

    #CATCH when the HLA does not exist yet
    #print(optitype_out_file)
    optitype_out_file = hla_files[0]
    if not os.path.exists(optitype_out_file):
        #print("WES WARNING: %s is not found!" % optitype_out_file)
        return ""

    f = open(optitype_out_file)
    hdr = f.readline().strip().split("\t") #ignore for now
    classI = f.readline().strip().split("\t")[1:7] #want first 6 cols
    #FOR classI alleles, prepend a HLA to each of them
    classI = ["HLA-%s" % a for a in classI if a]
    #print(classI)
    f.close()
    
    #check for xhla file
    classII = []
    if 'neoantigen_run_classII' in config and config['neoantigen_run_classII'] and len(hla_files) > 1:
        xhla_out_file = hla_files[1]
        
        #PARSE xhla json file...
        if os.path.exists(xhla_out_file):
            f = open(xhla_out_file)
            xhla_out = json.load(f)
            f.close()

            #build classII alleleles
            #ONLY add class II alleles--i.e. ones that start with "D"
            classII = [a for a in xhla_out['hla']['alleles'] if a.startswith("D")]
    if classII:
        classI.extend(classII)
    #NOTE: NOW classI has all hla alleles (including classII if opted for)
    hdr = ["A1", "A2", "B1", "B2", "C1", "C2"]
    hla = ["%s" % a for a in classI if a]
    if len(classI) > 6: #includes class II
        hdr.extend(["DP1","DP2","DQ1","DQ2","DR1","DR2"])
        hla = dict(zip(hdr, hla))
    else:
        hla = dict(zip(hdr, hla))
    #print(hla)
    return hla

###############################################################################
def getCNVInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    caller = config['somatic_caller']
    ret = []
    for run in config['runs']:
        #circos_plot = "wes_images/copynumber/%s.%s/circos.png" % (run, caller)
        sequenza_plot_1 = "wes_images/copynumber/%s/%s_cnv_genome.1.pdf" % (run,run)
        sequenza_plot_2 = "wes_images/copynumber/%s/%s_cnv_genome.2.pdf" % (run,run)
        sequenza_plot_3 = "wes_images/copynumber/%s/%s_cnv_genome.3.pdf" % (run,run)
        cnv_calls = "analysis/copynumber/%s/%s_cnvcalls.txt" % (run, run)
        tmp = {'name': run,
               #'circos_plot': circos_plot,
               'sequenza_plot_1': sequenza_plot_1,
               'sequenza_plot_2': sequenza_plot_2,
               'sequenza_plot_3': sequenza_plot_3,
               #Files
               'cnv_calls_file': (getFileName(cnv_calls), cnv_calls),
        }
        ret.append(tmp)
    #print(ret)
    return ret

def getPurityInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    #read in the purity file for each RUN
    ret = []
    for run in config['runs']:
        f_name = "analysis/purity/%s/%s.optimalpurityvalue.txt" % (run,run)
        if os.path.isfile(f_name):
            f = open(f_name)
            hdr = f.readline().strip().split("\t")
            vals = dict(zip(hdr, f.readline().strip().split("\t")))
            #print(vals)
            tmp = {'name': run, 'purity': vals['purity'],
                   'purity_file': (getFileName(f_name), f_name),
            }
            ret.append(tmp)
        else:
            print("WARNING: expected %s, but it does not exist" % f_name)
    return ret

def getHLAInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    ret = []
    for sample in config['samples']:
        optitype_fname = "analysis/optitype/%s/%s_result.tsv" % (sample,sample)
        xhla_fname = "analysis/xhla/%s/report-%s-hla.json" % (sample,sample)
        hla = parseHLA(config, [optitype_fname,xhla_fname])
        tmp = {'name': sample, 'HLA': hla,
               'optitype_file': (getFileName(optitype_fname), optitype_fname),
        }
        #Easy check for whether there were classII--check if it was set by
        #parseHLA
        if "DP1" in hla:
            tmp['xhla_file']= (getFileName(xhla_fname), xhla_fname)
        #print(tmp)
        ret.append(tmp)
        
    return ret

def getNeoantigenInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    ret = []
    for run in config['runs']:
        sub_path="wes_images/neoantigen/%s/" % run
        hla_plot = sub_path+"HLA_epitopes_fraction_plot.png"
        patient_plot = sub_path+"Patient_count_epitopes_plot.png"
        epitopes_plot = sub_path+"epitopes_affinity_plot.png"

        if 'neoantigen_run_classII' in config and config['neoantigen_run_classII']:
            pvacseq_file = "analysis/neoantigen/%s/combined/%s.all_epitopes.tsv" % (run, run)
        else:
            pvacseq_file = "analysis/neoantigen/%s/MHC_Class_I/%s.all_epitopes.tsv" % (run, run)
        tmp = {'name': run,
               'hla_plot': hla_plot,
               'patient_plot': patient_plot,
               'epitopes_plot':epitopes_plot,
                'pvacseq_file': (getFileName(pvacseq_file), pvacseq_file),
        }
        ret.append(tmp)
    return ret

def getClonalityInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    ret = []
    for run in config['runs']:
        density_plot = "wes_images/clonality/%s/%s_plot.density.png" % (run,run)
        scatter_plot = "wes_images/clonality/%s/%s_plot.scatter.png" % (run,run)
        coordinates_plot = "wes_images/clonality/%s/%s_plot.coordinates.png" % (run,run)
        table_file = "analysis/clonality/%s/%s_table.tsv" % (run,run)
        clonality = calcClonality(table_file)
        tmp = {'name': run,
               'density_plot': density_plot,
               'scatter_plot': scatter_plot,
               'coordinates_plot': coordinates_plot,
               'clonality': clonality,
               'pyclone_table_file': (getFileName(table_file), table_file),
        }
        ret.append(tmp)
    #print(ret)
    return ret

def main():
    usage = "USAGE: %prog -c [config.yaml file] -o [output html file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-c", "--config", help="config.yaml file")
    optparser.add_option("-o", "--output", help="output html file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.config or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #PARSE config.yaml
    config={}
    with open(options.config, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    config = getRuns(config)
    
    #ASSUMING it's being run as WES project level
    templateLoader = jinja2.FileSystemLoader(searchpath="cidc_wes/report")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template("wes_level2.html")

    #STANDARD nav bar list
    nav_list = [('wes_meta.html', 'WES_META'),
                ('wes_level1.html','WES_Level1'),
                ('wes_level2.html','WES_Level2')]

    #MODIFY this!
    pg_name = "WES_LEVEL_2"
    sidebar = [("cnv", "CNV", []),
               ("purity", "Purity", []),
               ('hla','HLA', []),
               ("neoantigen","Neoantigen", []),
               ("clonality","Clonality", []),
    ]
    
    wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                       'page_name': pg_name}

    #CNV
    cnv = getCNVInfo(config)
    wes_report_vals['cnv'] = cnv

    #PURITY
    purity = getPurityInfo(config)
    wes_report_vals['purity'] = purity

    #HLA
    hla = getHLAInfo(config)
    wes_report_vals['hla'] = hla

    #Neoantigen
    neoantigen = getNeoantigenInfo(config)
    wes_report_vals['neoantigen'] = neoantigen

    #Clonality
    clonality = getClonalityInfo(config)
    wes_report_vals['clonality'] = clonality

    template.stream(wes_report_vals).dump(options.output)  
        
if __name__ == '__main__':
    main()
