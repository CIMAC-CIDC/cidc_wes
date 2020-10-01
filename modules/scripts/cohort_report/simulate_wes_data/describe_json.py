#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser

def whatBasicType(elm):
    """IF elm is an int, float, or str, returns (True, <type>)"""
    if isinstance(elm, int):
        return(True, "int")
    elif isinstance(elm, float):
        return(True, "float")
    elif isinstance(elm, str):
        return(True, "string")
    else:
        return(False, None)
    
def describe(dict_a, indent=0):
    """Given an dict print the keys of the dictionary; if encounter an elem
    that is a dict, recursively call this fn to describe that dict
    """
    for k in dict_a.keys():
        elm = dict_a[k]
        space = " "*indent
        (basicType, ttype) = whatBasicType(elm)
        if basicType:
            print("%s%s:%s" % (space, k, ttype))
        elif isinstance(elm, list):
            #check first elm
            if len(elm) > 1:
                (isBasic, list_type) = whatBasicType(elm[0])
                if not isBasic:
                    list_type = "complex"
            print("%s%s:[%s]" % (space, k, list_type))
        elif isinstance(elm, dict):
            print("%s%s:" % (space, k))
            describe(elm, indent+4)
            
def main():
    usage = "USAGE: %prog -f json file to describe"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="json file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.file)
    tmp = json.load(f)
    f.close()

    describe(tmp)
    
if __name__=='__main__':
    main()
