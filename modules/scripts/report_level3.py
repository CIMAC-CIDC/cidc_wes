#!/usr/bin/env python3
"""Script to generate analysis/report/wes_level1.html"""

import os
import sys

import re
import subprocess
import yaml

from optparse import OptionParser
import jinja2
import pandas as pd

from report_level1 import getFileName
from report_level1 import getRuns
from clonality_calcClonality import calcClonality

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
    template = templateEnv.get_template("wes_level3.html")

    #STANDARD nav bar list
    nav_list = [('wes_meta.html', 'WES_META'),
                ('wes_level1.html','WES_Level1'),
                ('wes_level2.html','WES_Level2'),
                ('wes_level3.html','WES_Level3')]

    #MODIFY this!
    pg_name = "WES_LEVEL_3"
    sidebar = [("clonality", "Clonality", [])]
    
    wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                       'page_name': pg_name}

    #Clonality
    clonality = getClonalityInfo(config)
    wes_report_vals['clonality'] = clonality

    template.stream(wes_report_vals).dump(options.output)  
        
if __name__ == '__main__':
    main()
