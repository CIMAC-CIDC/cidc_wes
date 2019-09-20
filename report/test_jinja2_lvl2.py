#!/usr/bin/env python3
import jinja2

TEMPLATE_FILE = "wes_level2.html"
OUT_FILE = "output/wes_level2.html"


templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template(TEMPLATE_FILE)
#HERE is where we populate a dictionarity of associated values
nav_list = [('wes_level1.html','WES_Level1'),
            ('wes_level2.html','WES_Level2'),
            ('wes_level3.html','WES_Level3')]

sidebar = [("coverage", "Coverage"), ("cnv", "CNV"), ('purity','Purity'),
           ("hla","HLA"), ('neoantigen',"Neoantigen")]
pg_name = 'WES_LEVEL_2'

wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                   'page_name': pg_name}
    
template.stream(wes_report_vals).dump(OUT_FILE)  

