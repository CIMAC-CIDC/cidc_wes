
def getNormal(wildcards):
    tmp=[]

    print('helele')
    r = config['runs'][wildcards.run]
    print(r)

    #check that we have a valid pair
    if len(r) >=2:
        tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[0],r[0])]
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! Missing Normal sample for run: %s" % (wildcards.run)]
    print(tmp)
    return tmp

def getTumor(wildcards):
    tmp=[]
    
    r = config['runs'][wildcards.run]
    print(r)

    #check that we have a valid pair
    if len(r) >=2:
        tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[1],r[1])]
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! Missing Tumor sample for run: %s" % (wildcards.run)]
    print(tmp)
    return tmp

def targets(wildcards):
    """Generates the targets for this module"""
    #print(wildcards)
    ls = []
    for run in config["runs"].keys():
        ls.append("analysis/run_example/%s/%s.foobar.txt" % (run,run))
    return ls

rule example_all:
    input:
        targets
    
rule foo:
    """Dummy rule that does nothing ACTUALLY it breaks!"""
    input:
        norm = getNormal, #RETURNS SampleName under "Normal" in metasheet.csv
        tumor = getTumor #RETURNS SampleName under "Tumor" in metasheet.csv
    output:
        "analysis/run_example/{run}/{run}.foobar.txt"
    shell:
        "echo hello"
              
