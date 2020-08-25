#!/usr/bin/env python3
"""Script to generate wes json data"""

import os
import sys

import random
import json
import math

from optparse import OptionParser

_gc_content = [133008,96476,111922,142573,153722,156990,178073,191338,218133,246522,290547,355549,464525,613391,870301,1281550,1918346,2881743,4322293,6381663,9217796,12965431,17703536,23500602,30331714,38036948,46323412,54981608,63719368,72234510,80345995,87896316,94723078,100731539,105938290,110269883,113702244,115529374,115334607,113360183,109880066,104982262,99165281,93153305,87763393,83411311,80059055,77383561,74702671,71434843,67304414,62480720,57482129,52809488,48647914,44978082,41441569,37680657,33536258,29046522,24542771,20384107,16835666,13908082,11492178,9480842,7814208,6349664,5078099,4064637,3261633,2611527,2112841,1724151,1422238,1197298,1032749,909126,809173,706562,602758,507135,425573,349915,271605,197215,131726,85982,59302,38361,25087,16055,9563,6006,3424,1675,986,658,328,133,339]

_insert_size = [28,29,23,42,34,31,41,44,36,33,43,45,32,44,37,61,66,341,261,187,179,145,145,140,143,143,156,161,161,173,225,252,269,348,413,470,576,671,809,1051,1261,1446,1777,2143,2690,3119,3618,4462,5141,5908,6815,8165,9402,10636,12464,13982,15949,17995,20104,22622,25662,28666,31293,34757,38452,41957,46410,50419,55229,59539,64184,69603,74661,79941,86008,92608,98384,105001,111351,118200,124178,132625,139707,145642,154489,161833,169842,178042,185475,192404,201854,210156,217139,225742,235166,242862,250569,259750,267212,275239,284243,291419,300631,308028,317088,324158,332409,340457,348807,356359,362825,370942,377525,385271,391797,399538,406072,413332,419228,423996,432484,438167,444570,449129,454381,459595,465747,469727,474040,479066,484015,489060,492446,496954,500084,503061,508126,510972,513407,516923,519209,521066,523710,527294,528460,531407,530679,533921,534289,534535,536432,537664,538093,538956,537573,539214,538593,537581,536781,537614,535787,536263,536286,534760,533023,532367,531546,528278,526892,525708,523633,522639,520684,517393,516361,513482,510575,508138,506201,502698,500375,497057,493476,491685,489001,484624,481597,478612,475242,471521,469169,465211,461634,457718,453568,451364,446252,441872,438236,434659,430764,426735,424208,418227,415407,409711,405901,400672,398560,395345,390058,386398,380837,377880,373542,369334,365661,361547,356698,353305,349812,344559,339858,336475,331749,328100,324669,321663,314654,312548,307890,304593,299937,295397,293280,287233,283898,279976,275647,271883,267571,264776,260970,257322,253772,250829,246520,243532,240049,236478,232716,229370,225267,221783,219031,215897,212747,208959,206442,202759,199844,195675,193541,190529,187496,184791,181211,178705,176324,173221,170024,166924,165344,161751,160626,157808,154723,152071,149543,146937,144338,142684,140042,137311,135562,133499,130887,128398,126630,124495,122291,119972,117743,116232,113667,112703,110091,108101,106593,104685,102843,101832,99147,97465,96247,94544,92670,91228,88851,88114,86092,85038,83572,81644,80361,78271,77756,75876,74948,73632,72017,70328,69113,67979,66789,65981,64708,63548,62161,61341,59924,59025,57626,57003,55508,54934,53638,52941,52381,50922,50172,48776,47976,47113,46249,45615,44346,43822,42870,42187,41545,40782,39677,39240,38487,37558,36758,36119,35415,34839,34007,33424,32762,32329,31689,30983,30429,29925,29199,28274,27824,27515,26832,26559,26277,25590,24983,24340,23961,23505,23473,22883,21999,21590,21420,20872,20288,19981,19644,19500,18639,18159,18027,17751,17429,16902,16751,16254,16228,15653,15392,15073,14803,14471,14003,14043,13605,13107,12881,12721,12050,12231,11928,11741,11366,11360,11175,10750,10501,10115,10026,9827,9868,9414,9361,9258,8827,8644,8468,8290,8184,7986,7751,7655,7491,7461,7233,7043,6907,6826,6582,6319,6324,6060,5991,6005,5726,5673,5561,5505,5313,5211,5105,4945,4893,4772,4631,4633,4550,4335,4209,4040,4127,4003,3838,3770,3702,3634,3529,3573,3333,3348,3284,3190,3105,3083,3010,2932,2852,2828,2771,2679,2644,2582,2555,2486,2342,2375,2327,2297,2243,2043,1987,2029,1958,1945,1821,1855,1819,1784,1743,1647,1625,1602,1624,1598,1560,1523,1476,1428,1391,1414,1273,1316,1321,1225,1176,1212,1193,1191,1176,1091,1005,1062,1036,997,1023,914,916,955,918,926,903,839,859,845,864,793,785,751,734,682,707,718,674,668,699,673,642,646,632,578,603,590,587,536,550,530,506,537,466,496,466,442,487,439,436,464,465,440,404,409,378,391,376,391,347,375,360,357,321,300,336,345,292,288,289,271,263,275,282,247,248,248,266,265,266,267,227,213,236,212,241,231,224,200,197]

def makeAllele(base_str, range1, range2):
    """Given a base string, e.g. HLA-A, a range for the first set, and a range
    for the second set
    returns base_str*{randint(1, range1)}:{randint(1,range2))}
    NOTE: the random numbers are 0 padded
    """
    n = str(random.randint(1, range1)).zfill(2)
    m = str(random.randint(1, range2)).zfill(2)
    return "%s*%s:%s" % (base_str, n, m)

def generateSample():
    total_reads = random.randint(5*10**8, 6*10**8) #500 M and 600M
    less = random.randint(1*10**6, 10*10**6) #1M to 10 M
    mapped_reads = total_reads - less
    less = random.randint(1*10**6, 3*10**6) #1M to 3 M
    dedup_reads = mapped_reads - less

    sample_gc = list(map(lambda x: x + random.randint(int(-0.05*x),int(0.05*x)),
                         _gc_content))
    sample_is = list(map(lambda x: x + random.randint(int(-0.05*x),int(0.05*x)),
                         _insert_size))
    alignment = {'total_reads': total_reads,
                 'mapped_reads': mapped_reads,
                 'dedup_reads': dedup_reads,
                 'gc_content': sample_gc,
                 #TO ADD
                 'mean_quality_score': random.gauss(36.0, 5),
                 'insert_size': sample_is}
    mean_depth= random.gauss(200.0, 50)
    median_depth= random.randint(int(mean_depth - 100), int(mean_depth))
    q1 = int(median_depth - (2.0/3.0)*median_depth)
    q3 = median_depth*2
    percent_bases = random.randint(70,100)
    
    coverage = {'total_reads': alignment['total_reads'],
                'mean_depth': mean_depth,
                'median_depth': median_depth,
                'q1_depth': q1,
                'q3_depth': q3,
                'percent_bases_gt_50': percent_bases}
    hla = {"A-1": makeAllele("HLA-A", 80, 10),
           "A-2": makeAllele("HLA-A", 80, 10),
           "B-1": makeAllele("HLA-B", 83, 10),
           "B-2": makeAllele("HLA-B", 83, 10),
           "C-1": makeAllele("HLA-C", 18, 10),
           "C-2": makeAllele("HLA-C", 18, 10),
           "DPB1-1": makeAllele("DPB1", 720, 1),
           "DPB1-2": makeAllele("DPB1", 720, 1),
           "DQB1-1": makeAllele("DQB1", 6, 99),
           "DQB1-2": makeAllele("DQB1", 6, 99),
           "DRB1-1": makeAllele("DRB1", 16, 48),
           "DRB1-2": makeAllele("DRB1", 16, 48),
           }
    return {'alignment': alignment, 'coverage': coverage, 'hla': hla}

def generateCopyNumber():
    "returns a dictionary of simulated cnv fields"
    clonality = round(random.random(), 3)
    purity = round(random.random(), 3)
    ploidy = random.gauss(2.5, 0.5) #mean 2.5, stdev 0.5
    ploidy = round(ploidy, 4)
    dipLogR = round(random.gauss(0, 1), 4)

    return {'clonality': clonality, 'purity': purity, 'ploidy': ploidy,
            'dipLogR': dipLogR,
            #HARD-CODED FILES--Stil unclear how to handle them
            'cnv_file': "analysis/clonality/CTTTP07T1.00/CTTTP07T1.00_pyclone.tsv",
            'cnv_plot_file': "analysis/report2/copy_number/01_copynumber_plot.png"}

def generateRunName():
    center = random.choice(['broad', 'mda', 'mocha'])
    sample = random.randint(0,1000)
    return "%s_s%s" % (center, sample)

def generateSomatic():
    """returns a dictionary of simulated somatic fields
    TODO: filtered_vcf_file, filtered_maf_file, tmb, functional_summary,
    trinucleotide_matrix, transition_matrix
    """
    snp = round(random.gauss(1000, 150)) #mean 1k, stdev 150
    insert = random.randint(5, 100)
    delete = random.randint(10, 100)
    total = snp + insert + delete

    mut_sum = {'snp': snp, 'insertion': insert,
               'deletion': delete, 'total': total}

    tumor = round(random.gauss(45000, 1500)) #mean 45k, stdev 1.5k
    normal = round(random.gauss(45000, 1500)) #mean 45k, stdev 1.5k
    overlap = round(random.random(), 3)
    common = round(tumor*overlap)

    tmb = {'tumor': tumor, 'normal': normal, 'common': common, 'overlap': overlap}

    ttypes = ['missense', 'nonsense', 'silent']
    func_summary = {}
    for t in ttypes:
        snp = round(random.gauss(100, 25))
        insert = random.randint(5, 50)
        delete = random.randint(10, 100)
        foo = {'snp': snp, 'insert': insert, 'delete': delete, 'total': snp + insert + delete}
        func_summary[t] = foo
        
    #Transition matrix
    ACGT = ["A","C","G","T"]
    trans_mat = {}
    for k in ACGT:
        tmp = dict([(l,random.randint(50,350)) for l in ACGT])
        tmp[k] = 0
        trans_mat[k] = tmp
    #print(trans_mat)
    return {'mutation_summary': mut_sum, 'transition_matrix': trans_mat, 'tmb':tmb, 'functional_summary': func_summary}

def generateMeta(run_name):
    """returns a dictionary of simulated meta fields- e.g. age, sex, etc
    NOTE: like viper, this can be any abitrary variable that can be added
    to the cohort metasheet
    """
    #parse the runname to get the center
    center = run_name.split("_")[0]

    age = random.randint(20,75)
    sex = random.choice(['M', 'F'])
    treatment = random.choice(['A', 'B'])
    response = random.choice(['Responder', 'Non-Responder'])
    return {'center': center, 'sex': sex, 'age': age, 'treatment':treatment,
            'response': response}
    
def main():
    usage = "USAGE: %prog -s [seed] -o [output json file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-s", "--seed", help="random seed")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)


    if options.seed:
        random.seed(options.seed)

    run_name = generateRunName()
    meta = generateMeta(run_name)
    tumor = generateSample()
    tumor['id'] = "%s_T" % run_name
    normal = generateSample()
    normal['id'] = "%s_N" % run_name

    cnv = generateCopyNumber()
    #print(cnv)
    somatic = generateSomatic()
    #print(somatic)

    run = {'id': run_name,
           'meta': meta,
           'tumor': tumor,
           'normal': normal,
           'copy_number': cnv,
           'somatic': somatic,
    }
    
    if options.output:
        out = open(options.output,'w')
    else:
        #write to the run_name.json
        out = open("%s.json" % run_name, 'w')
        
    out.write(json.dumps(run))
    out.close()
    
if __name__ == '__main__':
    main()
