#!/usr/bin/env python3
"""Len Taing 2020 (TGBTG)
Generalized script to generate analysis/report/report.html
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

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

def is_image_file(s):
    """Checks to see if the string starts with 'img:'"""
    return s.startswith('img:')

def chop_image_file(s):
    return s[4:]

def captionHelper(extension, path, basename):
    """Given an extenstion like 'caption.txt' or 'sub_caption.txt', 
    a base path and a base file name, tries to look to see if there is 
    a caption/subcaption--returns the contents of the file if it exists
    otherwise None"""
    caption = None
    if os.path.exists(os.path.join(path, "%s_%s" % (basename,extension))):
        f = open(os.path.join(path, "%s_%s" % (basename,extension)))
        caption = f.read().strip()
        f.close()
    return caption

#NOTE: this can easily handle csv files too!
def buildTable(tsv_file, jinjaEnv, separator):
    """Given a tsv file, and a section--converts the tsv file to a table
    assumes the first line is the hdr"""
    jinjaEnv.tests['imagefile'] = is_image_file
    jinjaEnv.filters['chopimagefile'] = chop_image_file
    template = jinjaEnv.get_template("table.html")
    #Title is the filename prettyprinted
    fname = ".".join(tsv_file.split("/")[-1].split(".")[:-1])
    path = "/".join(tsv_file.split("/")[:-1]) #drop the file
    title = prettyprint(fname, True)

    vals = {'title':title }

    #Check for a caption
    caption = captionHelper('caption.txt', path, fname)
    if caption:
        vals['caption'] = caption
    #check for subcaption
    sub_caption = captionHelper('subcaption.txt', path, fname)
    if sub_caption:
        vals['sub_caption'] = sub_caption
        
    table = []
    f = open(tsv_file)
    hdr = f.readline().strip().split(separator)
    
    for l in f:
        tmp = l.strip().split(separator)
        table.append(tmp)
    vals['header'] = hdr
    vals['table'] = table

    #print(vals)
    return template.render(vals)

def buildPlot(png_file, jinjaEnv):
    """Given a png file displays the plot..simple!"""
    template = jinjaEnv.get_template("plot.html")
    fname = ".".join(png_file.split("/")[-1].split(".")[:-1])
    path = "/".join(png_file.split("/")[:-1]) #drop the file
    title = prettyprint(fname, True)

    #make the png file path REALITIVE to the report.html file!
    png_file_relative = "/".join(png_file.split("/")[2:])
    vals = {'title':title,
            'png_file': png_file_relative
    }
    #Check for a caption
    caption = captionHelper('caption.txt', path, fname)
    if caption:
        vals['caption'] = caption
    #check for subcaption
    sub_caption = captionHelper('subcaption.txt', path, fname)
    if sub_caption:
        vals['sub_caption'] = sub_caption
        
    #print(vals)
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
    first_section = ""
    #ONLY sections so far--no further recursion
    
    #for sect in os.listdir(options.dir):
    for (i, sect) in enumerate(sections):
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
        if i == 0: #first element is shown
            first_section = sect
            tmp = """<div id="%s" class="container wes_container">\n""" % sect
        else:
            tmp = """<div id="%s" class="container wes_container" style="display:none">\n""" % sect
            
        for ffile in ordering:
            filepath = os.path.join(path, ffile)
            if os.path.isfile(filepath): #SKIP directories
                #CHECK for file existance
                if not os.path.exists(filepath):
                    print("WARNING: report.py- file %s is not found. SKIPPED for rendering" % filepath)
                    continue
                if ffile.endswith(".tsv"): #MAKE a table
                    tmp += buildTable(filepath, templateEnv, '\t')
                elif ffile.endswith(".csv"):
                    tmp += buildTable(filepath, templateEnv, ',')
                elif ffile.endswith(".png"): #Make a plot
                    tmp += buildPlot(filepath, templateEnv)
        #END container
        tmp += "\n</div>"
        wes_panels[sect] = tmp

    wes_sections = [(s, prettyprint(s)) for s in sections]
    wes_report['sections'] = wes_sections
    wes_report['panels'] = wes_panels
    wes_report['first_section'] = first_section
    template.stream(wes_report).dump(options.output)  

if __name__ == '__main__':
    main()
