#!/usr/bin/env python3
import jinja2

TEMPLATE_FILE = "wes_level1.html"
OUT_FILE = "output/wes_level1.html"


templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(TEMPLATE_FILE)
#HERE is where we populate a dictionarity of associated values
nav_list = [('wes_level1.html','WES_Level1'),
            ('wes_level2.html','WES_Level2'),
            ('wes_level3.html','WES_Level3')]
sidebar = [("meta", "META", []),
           ("alignment", "Alignment", ['Mapping_Stats', 'GC_Bias_Plots','Quality_Score','Quality_by_Cycle']),
           ('coverage','Coverage', []),
           ("somatic","Somatic", []),
           ('germline',"Germline", [])]
pg_name = 'WES_LEVEL_1'

meta = {'wes_version' : "v1.1 (commit: d8c124c)",
        'ref_version' : "ver1.0 (build date: 20190911)",
        "assembly_version": "GDC hg38",
        "sentieon_version": "201808.05",
        "somatic_caller": "tnscope (sentieon)",
        "neoantigen_callers": "MHCflurry NetMHCcons MHCnuggetsII",
        "epitope_lengths": "8,9,10,11",
        "snakemake_version": "5.4.5"}
wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                   'page_name': pg_name}
#INSERT meta values into the dictionary!
for k in meta.keys():
    wes_report_vals['meta_%s' % k] = meta[k]
    
template.stream(wes_report_vals).dump(OUT_FILE)  

