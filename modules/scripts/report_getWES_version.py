#!/usr/bin/env python
"""Script to get the WES git commit string"""

import os
import sys
import subprocess
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -o [output file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-o", "--out", help="output file to store wes version", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.out:
        optparser.print_help()
        sys.exit(-1)

    #CHANGE to cidc_wes, but first store this
    wd = os.getcwd()
    os.chdir("cidc_wes")
    #GET wes commit string--first six letters of this cmd
    cmd = "git show --oneline -s".split(" ")
    output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    wes_commit = output[:6].decode("utf-8") 

    #GET wes current tag
    cmd = "git describe --tags".split(" ")
    output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    wes_tag = output.decode("utf-8").split("-")[0]
    wes_version = "WES %s (commit: %s)" % (wes_tag, wes_commit)
    #CHANGE back to CWD
    os.chdir(wd)


    #print(wes_version)
    #print(wes_tag)
    #print(output.decode("utf-8"))

    out = open(options.out, "w")
    out.write("%s\n" % wes_version)
    
if __name__=='__main__':
    main()
