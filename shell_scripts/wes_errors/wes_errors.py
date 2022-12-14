#!/usr/bin/env python

#Jacob Geisberg 8/2022
import pydot
import pandas as pd
import requests
import json
import os
import sys
import argparse
import yaml
from pathlib import Path
import subprocess
#from collections import deque

code_prompt= ("error codes:\n"
'00: Unknown error\n'
'01: Data file corrupt or unreadable\n'
'02: PE reads incorrectly paired\n'
'03: Poor read quality\n'
'04: Too few reads\n'
'05: Duplicate normal sample -- not used for analysis\n'
'11: External software bug/issue {software: XXX}\n'
'12: WES pipeline software bug/issue {module: XXX}\n'
'13: Upstream module bug/issue {module: XXX}\n'
'21: Out of disk space error\n'
'22: Out of memory error\n'
'23: Unexpected interruption of service\n'
'30: Did not finish computation after 6 hours\n')

codes={
    '00': 'Unknown error',
    '01': 'Data file corrupt or unreadable',
    '02': 'PE reads incorrectly paired',
    '03': 'Poor read quality',
    '04': 'Too few reads',
    '05': 'Duplicate normal sample -- not used for analysis',
    '11': 'External software bug/issue {software: XXX}',
    '12': 'WES pipeline software bug/issue {module: XXX}',
    '13': 'Upstream module bug/issue {module: XXX}',
    '21': 'Out of disk space error',
    '22': 'Out of memory error',
    '23': 'Unexpected interruption of service',
    '30': 'Did not finish computation after 6 hours'
    }

def get_downstream(node, nodes):
    '''takes a job ID and identifies downstream job IDs, including self, effectively a BFS'''
    downstream = [i for i in nodes[node]]
    downstream.append(node) # WE NEED TO HANDLE BROKEN RULE TOO!
    queue = [i for i in nodes[node]]
    while queue:
        n = queue.pop()
        downstream.append(n)# we will get unique values later
        if n in nodes: # check if node has children since some rules are endpoints
            children=nodes[n]
            for child in children:
                if child not in queue:
                    queue.append(child)
    return list(set(downstream))




def get_node_files(node, graph, files):
    '''takes a job number, job graph, and job to file dataframe and returns a list
    of the files prodcued by the job. '''

    n = graph.get_node(node)[0]
    label=n.get("label").split("\\n")[0].strip('\"')
    f=files.loc[files['rule'] == label, ['output_file']] #match label of job to rule in files
    f=f.output_file.to_list()
    return f

def get_ingested_files(id, TO, tumor, normal=None):
    '''Pulls all possible ingested files from dropbox and returns the names of all
    files needed for the run based on runID and tunor and normal names. This list
    includes files that have already been sucessfully run.'''

    print("getting ingetsed files...")

    #PULLING FILES FROM GITHUB
    url = 'https://raw.githubusercontent.com/CIMAC-CIDC/cidc-ngs-pipeline-api/master/cidc_ngs_pipeline_api/wes/wes_output_API.json'
    files = requests.get(url)
    files = files.json()
    wildcards = {"run id": id,
                 "tumor cimac id": tumor,
                 "normal cimac id": normal}

    # COMPILING RELEVANT FILE NAMES BASED ON SUPPLIED NAMES AND DATA
    lst=[]
    for wildcard in files.keys(): # API is split by run_id, tumor and normal files
        if not ((wildcard == "normal cimac id") and TO == True): # dont make normal sample files for TO runs. API may need fixing
            for file in files[wildcard]:
                path = file['file_path_template']
                file_TO = file['tumor_only_assay']
                optional = file['optional']
                exclude = TO and (not file_TO)#ensures that no normal files are included for TO samples
                path = path.replace('{'+ wildcard +'}', wildcards[wildcard])
                if (not optional) and (not exclude):
                    lst.append(path)
    return lst

def get_absolute_path(file, folder):
    #DETERMINE ABSOLUTE FILE PATH
    if folder[-1] == "/":
        path = folder + file
    else:
        path = folder + "/" + file
    return path


def file_writer(path):
    '''Writes text to a given path. Will create directories as needed.
    DOES NOT OVERWITE EXISTING FILES!
    '''

    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not os.path.exists(path):
        Path(path).touch()
        print('wrote: %s' %(path))
        # # for testing purposes only
        # with open(path, "w") as f:
        #     f.write("I am Jason Bourne")




def main():
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, # keep newlines
    description="Finds all affected files from a given error, generates blank files that are needed for ingestion, and allows easy construction of error.yaml",
    epilog="")
    group = parser.add_mutually_exclusive_group(required=True)
    #commenting these inputs until they are implemented
    #group.add_argument("-r", "--rule", help="name of the rule that errored out")
    #group.add_argument("-f", "--file", help="name of the output file that errored out")
    group.add_argument("-j", "--job", help="job id number that errored out")
    parser.add_argument("--folder", default="/mnt/ssd/wes/", help="location of WES results. default: /mnt/ssd/wes/")# should it be /mnt/ssd/wes/analysis?
    parser.add_argument("-d", "--dag", default="sample_dag.dot", help="file containing a complete dag of the run. can be generated with: snakemake -s cidc_wes/wes.snakefile --forceall -n --dag > sample_dag.dot")
    parser.add_argument("-m", "--metasheet", default="metasheet.csv", help="name of the metasheet used. default: metasheet.csv")
    args = parser.parse_args()
    argument_dict = vars(args) #needed?


    #GET RUN, TUMOR AND NORMAL NAMES
    samples=pd.read_csv(args.folder + "/" + args.metasheet, header=15, keep_default_na=False).loc[0]
    run_name = samples["RunName"]
    tumor = samples["Tumor"]
    normal = samples["Normal"]

    if normal == '':
        TO = True
        dag = 'tumor_only_dag.dot'
        summary = 'tumor_only_summary.tsv'
    else:
        TO = False
        dag = 'tumor_normal_dag.dot'
        summary = 'tumor_normal_summary.tsv'


    #CREATE JOB GRAPH
    graph = pydot.graph_from_dot_file(dag)
    graph=graph[0]
    nodes={}
    for edge in graph.get_edge_list():
        source=edge.get_source()
        destination=edge.get_destination()
        if source not in nodes:
            nodes[source]=[destination]
        else:
            nodes[source].append(destination)

    # READ FILE TO RULE TABLE
    files=pd.read_csv(summary, header=1, delimiter="\t")


    #GET DOWNSTREAM MODULES AND FILES (RIGHT NOW ONLY JOBID IS SUPPORTED)
    #SOME PSEUDO-CODE HERE FOR ALTERNATIVE INPUTS
    #if rule loop over nodes in dag and return job of rule in # qeustion
    # if file, parse file, get rule, get job
    #print(args.job)
    downstream_jobs=get_downstream(args.job, nodes)
    downstream_files = []
    for job in downstream_jobs:
        downstream_files += get_node_files(job, graph, files)

    #CONVERT FILES INTO PROPER NAMES BASED ON METASHEET
    for i in range(len(downstream_files)):
        downstream_files[i] = downstream_files[i].replace("RUN", run_name)
        downstream_files[i] = downstream_files[i].replace("TUMOR", tumor)
        downstream_files[i] = downstream_files[i].replace("NORMAL", normal)
        downstream_files[i]
    print()




    #GET INGESTED FILES FROM CIDC API AND INTERSECT WITH DOWNSTREAM FILES
    ingested_files=get_ingested_files(run_name, TO, tumor, normal)
    #we use a set comprehension to create a sorted set of the intersection of downstream and ingested files
    affected_files=sorted({x for x in downstream_files if x in ingested_files})
    print('The following files are affected by the error and are nessisary for ingestion:')
    for x in affected_files:
        print(x)
    print()

    #CHECKS ALL INGESTED FILES FOR ANY MISSING ONES, AND PROMPTS USER IF WE FIND THEM
    unacconted = [f for f in ingested_files if not os.path.exists(get_absolute_path(f, args.folder))]
    unacconted = [f for f in unacconted if f not in affected_files]
    if len(unacconted) > 0:
        print("The following required files are missing but are not downstream of the selected error:")
        for f in unacconted:
            print(f)
        proceed = False
        while proceed not in ['y', 'n']:
            proceed = input("Would you like to add these files to the affected list and proceed? (y/n): ")
            if proceed == "y":
                affected_files.extend(unacconted)
            elif proceed == 'n':
                print("Aborting! Nothing was written.")
                sys.exit(-1)



    #HANDLE EACH FILE
    yaml_dict = {"errors":{}}
    previous=""
    for file in affected_files:
        path = get_absolute_path(file, args.folder)
        print("CURRENT FILE:", path)

        #ADD FILE STATUS TO HELP THE USER
        if not os.path.exists(path):
            print("FILE STATUS: missing")
            file_writer(path) # we create missing files for ingestion

        elif os.path.getsize(path) == 0:
            print("FILE STATUS: empty")

        #add condition here to handle unhealthy looking file

        else:
            print("FILE STATUS: uncertain")
            #if file is not empty, we add a preview
            print("FILE PREVIEW:")
            subprocess.run("head -n 5 %s" % (path), shell=True)
            #print()
            # if os.path.exists(path):
            #     print("file_preview:")
            #     subprocess.run("head -n 5 %s" % (path), shell=True)
            #     print()
            print("END FILE PREVIEW")


        #GETTING ERROR CODE FROM USER
        valid=False
        while not valid:
            error_code = input("Enter a valid error code. Press 'e' to view error codes, 'f' to view file_path, and 'ENTER' to add no error: ")
            if error_code == "":
                valid = True
                print("no error for %s" % path)
            elif error_code == "e":
                print(code_prompt)
            elif error_code == "f":
                print(file)
            elif error_code not in codes:
                print('Error code invalid. Please try again!')
            else:
                valid = True
                if error_code in ["11","12","13"]:
                    if error_code == "11":
                        prompt = "Please enter errored software for the file: "
                    elif error_code == '12':
                        prompt = "Please enter errored module for the file: "
                    else:
                        prompt = "Please enter the upstream module error for the file: "
                        # possibly include autocomplete or memory function here but there are downsides
                    txt =  input(prompt)

                #ADDING ERROR MESSAGE TO OUTPUT DICTIONARY
                message = "ERROR%s: %s" % (error_code, codes[error_code])
                if error_code in ["11","12","13"]:
                    message = message.replace("XXX", txt)
                yaml_dict['errors'][file] = [{"error": message}]
                print(file + ':', yaml_dict['errors'][file]) # testing code


        #TAKING COMMENTS FROM THE USER
        valid = False
        while not valid:
            comment = input("Enter a comment, press 'p' to repeat previous comment, or press 'Enter' to skip: ")
            if comment == '':
                valid =True
            else:
                if comment == "p":
                    comment = previous
                response=input("You are about to add '%s' as a comment to %s. Is this ok (y/n): " % (comment, file))
                if response == "y":
                    valid = True
                    previous = comment

                    #ADDING COMMENT TO OUTPUT DICTIONARY
                    if file in yaml_dict['errors']:
                        yaml_dict['errors'][file].append({"comments": comment})
                    else:
                        yaml_dict['errors'][file]=[{"comments": comment}]
        print()# for spacing


    #WRITE YAML
    p = args.folder + '/' + "analysis/%s_error.yaml" % (run_name)
    with open(p, "w") as outfile:
        yaml.dump(yaml_dict, outfile, width=float("inf"), sort_keys=True)

    print("successfully wrote yaml to: %s" % (p))


if __name__=='__main__':
    main()
