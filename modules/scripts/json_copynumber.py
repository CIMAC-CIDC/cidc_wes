#!/usr/bin/env python
"""
Len Taing 2022 (TGBTG)
Given cnvkit, cnvkit_ehanced, sequenza, facets, consensus, consensus_merged_gain, and consensus_merged_loss files, encodes theses in base64 and stashes them
"""

import os
import sys
import subprocess
import json
import base64

from optparse import OptionParser

def base64encode(in_file):
    """encodes the contends of a file and returns a b64 string"""
    #ENCODE the input and results files
    f = open(in_file)
    s = f.read()
    s_byte = s.encode('utf-8')
    #NOTE: .decode is needed to convert the bytes back to a string
    s_b64 = base64.b64encode(s_byte).decode('utf-8')
    f.close()
    
    return s_b64
    

def main():
    usage = "USAGE: %prog -r run_name -f file -o output_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("--cnvkit", help="raw cnvkit results file to stash", default=None)
    optparser.add_option("--cnvkit_enhanced", help="raw cnvkit enhanced results file to stash", default=None)
    optparser.add_option("--sequenza", help="raw sequenza enhanced results file to stash", default=None)
    optparser.add_option("--facets", help="raw facets enhanced results file to stash", default=None)
    optparser.add_option("--consensus", help="consensus cnv bed file to stash", default=None)
    optparser.add_option("--consensus_gain", help="consensus merged gains cnv bed file to stash", default=None)
    optparser.add_option("--consensus_loss", help="consensus merged loss cnv bed file to stash", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.cnvkit or not options.cnvkit_enhanced or not options.sequenza or not options.facets or not options.consensus or not options.consensus_gain or not options.consensus_loss or not options.output:
        optparser.print_help()
        sys.exit(-1)

    cnvkit = base64encode(options.cnvkit)
    cnvkit_enhanced = base64encode(options.cnvkit_enhanced)
    sequenza = base64encode(options.sequenza)
    facets = base64encode(options.facets)
    consensus = base64encode(options.consensus)
    consensus_gain = base64encode(options.consensus_gain)
    consensus_loss = base64encode(options.consensus_loss)
    
    js_out = {'id': options.run, 'copy_number': {'cnv': {'cnvkit': cnvkit, 'cnvkit_enhanced': cnvkit_enhanced, 'sequenza': sequenza, 'facets': facets, 'consensus': consensus, 'consensus_gain': consensus_gain, 'consensus_loss': consensus_loss}}}
    
    out = open(options.output, 'w')
    out.write(json.dumps(js_out))
    out.close()
    
if __name__=='__main__':
    main()
