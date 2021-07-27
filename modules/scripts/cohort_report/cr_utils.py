"""Len Taing (TGBTG)
Script of helper fns"""

import os
import sys
import json

def processSampleJson(json_fpath, subpart, attrs):
    """Reads in the WES json file, and will try to extract the attributes in 
    the attrs list from the subpart.
    NOTE: subparts must be in the set {'alignment','coverage', 'hla'}
    NOTE: normal samples information may be missing!

    returns a dictionary of values {'id': run_id, 'tumor': {tumor sample values}, 'normal': {normal sample values}}
    """
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()
    run = {'id': tmp['id']}

    #PROCESS tumor- based on old getSampleInfo fn
    sample = tmp['tumor']
    sub = sample[subpart]
    tumor = {'id': sample['id']}
    for a in attrs:
        tumor[a] = sub.get(a)
    run['tumor']= tumor

    #PROCESS normal- based on old getSampleInfo fn
    if 'normal' in tmp: 
        sample = tmp['normal']
        sub = sample[subpart]
        normal = {'id': sample['id']}
        for a in attrs:
            normal[a] = sub.get(a)
        run['normal']= normal
    return run

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

def millify(n):
    """Given a large int n, returns a string representation of n in human-
    readable form
    ref: https://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
    """
    millnames = ['',' K',' M',' B',' T']

    n = float(n)
    millidx = max(0,min(len(millnames)-1,
                    int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.1f}{}'.format(n / 10**(3 * millidx), millnames[millidx])
