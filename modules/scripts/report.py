#!/usr/bin/env python3
"""Generalized script to generate analysis/report/report.html
based on the contents of analysis/report sub-dirs: each sub-dir is a 'section'
"""

import os
import sys

import re
import subprocess
import yaml

from optparse import OptionParser
import jinja2
import pandas as pd

def parseYaml(yaml_file):
    """Parses a yaml file and returns a dictionary"""
    ret = {}
    with open(yaml_file, 'r') as stream:
        try:
            ret = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return ret

def prettyprint(s):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    return s.title()

def buildTable(tsv_file, section, jinjaEnv):
    """Given a tsv file, and a section--converts the tsv file to a table
    assumes the first line is the hdr"""
    template = jinjaEnv.get_template("table.html")
    #Title is the filename prettyprinted
    title = prettyprint(tsv_file.split("/")[-1].split(".")[0])

    vals = {'title':title,
            'section':section,
    }
    table = []
    f = open(tsv_file)
    hdr = f.readline().strip().split("\t")
    
    for l in f:
        tmp = l.strip().split("\t")
        table.append(tmp)
    vals['header'] = hdr
    vals['table'] = table
    return template.render(vals)

def main():
    usage = "USAGE: %prog -o [output html file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-d", "--dir", help="report directory path")
    optparser.add_option("-s", "--sections", help="sections list in order of appearance")
    optparser.add_option("-o", "--output", help="output html file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.dir or not options.sections or not options.output:
        optparser.print_help()
        sys.exit(-1)
    
    #ASSUMING it's being run as WES project level
    templateLoader = jinja2.FileSystemLoader(searchpath="cidc_wes/report2")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template("index.html")

    #Build up this dictionary
    wes_report = {}
    
    #infer sections from the analysis/report dir structure
    #sections = os.listdir(options.dir)
    sections = options.sections.split(",")
    wes_panels = {}

    #ONLY sections so far--no further recursion
    
    for sect in os.listdir(options.dir):
        if sect == "static": #SKIP static content if it's there
            continue
        #Check for {section}.yaml file for overrides
        path = os.path.join(options.dir, sect)
        overrides = {}
        #check for .yaml and .yml
        if os.path.exists(os.path.join(path, "%s.yaml" % sect)):
            overrides = parseYaml(os.path.join(path, "%s.yaml" % sect))
        elif os.path.exists(os.path.join(path, "%s.yml" % sect)):
            overrides = parseYaml(os.path.join(path, "%s.yml" % sect))

        #print(overrides)
        #Determine panel ordering--either from user override or default
        #which is sorted file list
        ordering = overrides.get("order",None) if overrides else None
        if not ordering:
            ordering = sorted(os.listdir(path))

        #Build container
        tmp = """<div id="%s" class="container wes_container">\n""" % sect
        for ffile in ordering:
            filepath = os.path.join(path, ffile)
            #CHECK for file existance
            if not os.path.exists(filepath):
                print("WARNING: report.py- file %s is not found. SKIPPED for rendering" % filepath)
                continue
            if ffile.endswith(".tsv"): #MAKE a table
                tmp += buildTable(filepath, sect, templateEnv)
        #END container
        tmp += "\n</div>"
        wes_panels[sect] = tmp

    wes_sections = [(s, prettyprint(s)) for s in sections]
    wes_report['sections'] = wes_sections
    wes_report['panels'] = wes_panels
    template.stream(wes_report).dump(options.output)  

if __name__ == '__main__':
    main()
