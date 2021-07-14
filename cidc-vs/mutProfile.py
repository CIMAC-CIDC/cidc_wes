#!/usr/bin/env python3
''' Mutation profile

usage: mutProfile.py [-h] -o OUTPUT -r REF -m MAF [-c CANCER] [-n NAME]

Tools to download public genomic data

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT (Plot)
                        Output Path (default: None)
  -j JSON OUTPUT --json_out output of trinucleotide context counts (json)
  -r REF, --ref REF     Reference Fasta file (default: None)
  -m MAF, --maf MAF     MAF File (default: None)
  -c CANCER, --cancer CANCER 
                        Reference cancer mutation frequency file (default: None )
  -n NAME, --name NAME  Sample Name (default: sample)
'''
import matplotlib.pylab as plt
import pandas as pd
from collections import OrderedDict
from io import StringIO
import os
import sys
import numpy as np
import argparse
import pybedtools
import json

def fetchMut(maf, only_protein=True, dna_alt_col="HGVSc", mut_type='snv'):
    """ Fetch the snv columns in maf file

    Parameters
    ----------
    maf : str
        /path/to/maf_file
    only_protein : bool, optiona
        Whether only extract snv at protein coding region (the default is True)
    dna_alt_col : str
        name of nucleotide change in the maf file (the default is "cDNA_Change")

    Returns
    -------
    pd.DataFrame
        A data frame only contain following information:
            Gene symbol
            mutated position (+/-1 extened, for 96 mutation profile)
    """

    df = pd.read_table(maf, comment='#')
    keep_col = ["Hugo_Symbol", "Chromosome",
                "Start_Position", "End_Position", dna_alt_col]

    if mut_type == 'snv':
        df = df.loc[(df["Variant_Type"] == "SNP") &
                    (~df[dna_alt_col].isna()), :]
    elif mut_type == 'indel':
        df = df.loc[(df["Variant_Type"].isin(['INS', 'DEL'])) &
                    (~df[dna_alt_col].isna()), :]
    else:
        raise KeyError(
            '{} is not valid, Only accept: indel, snv.'.format(mut_type))
    # extend mutated region
    df['End_Position'] = df['End_Position'] + 1
    df['Start_Position'] = df['Start_Position']-2

    if only_protein == True:
        return df.loc[~df["Variant_Classification"].isin(["3'Flank",
                                                          "3'UTR",
                                                          "5'Flank",
                                                          "5'UTR",
                                                          "Intron",
                                                          "RNA",
                                                          "IGR",
                                                          "Splice_Region",
                                                          "Splice_Site",
                                                          "Translation_Start_Site"]), keep_col]
    else:
        return df[keep_col]

def genTrinucleotideMtrx(maf, ref_fasta, dna_alt_col="HGVSc"):
    ''' Get somatic muation frequency based on maf file. 
    1. Extract sequencce from reference fasta file on corresponding mutation position.
    2. Count mutation frequency
  
    Parameters
    ----------
    maf : str
      Path to maf file
    ref_fasta : str
      Path to reference file 
    dna_alt_col : str, optional
      Column name of mutation indicator(e.g c.708G>C)(the default is "HGVSc", which [default_description])
  
    Returns
    -------
    pandas.DatFrame 
      Mutation frequency data frame with 
      - columns names: ['Neighbor', 'Alt', "Seq"] 
        - Neigbor: Nucleotide nearby snv. e.g A_C 
        - Alt: SNV
        - Seq: frequency of mutation grouped by Alt and Neighbor
      - index: SNV ["T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G"] 

    '''

    conv = dict(zip(("T>C", "A>G", "C>T", "G>A", "T>A", "A>T","T>G", "A>C", "C>A", "G>T", "C>G", "G>C"),
                    ("A>G", "A>G", "C>T", "C>T", "A>T", "A>T","A>C", "A>C", "C>A", "C>A", "C>G", "C>G")))

    snv_maf = fetchMut(maf, dna_alt_col=dna_alt_col).iloc[:, 1:]
    snv_intervel = pybedtools.BedTool.from_dataframe(snv_maf)
    ref_fasta = pybedtools.example_filename(ref_fasta)
    snv_obj = snv_intervel.sequence(fi=ref_fasta, tab=True, name=True)
    snv_seq = StringIO(open(snv_obj.seqfn).read())
    
    df = pd.read_table(snv_seq, names=[dna_alt_col, 'Seq'])
    df['Neighbor'] = df.Seq.map(lambda x: x[0]+"_"+x[2])
    df['Alt'] = df[dna_alt_col].map(lambda x: x[-3:]).map(conv)
    return df.groupby(['Neighbor', 'Alt']).count().reset_index()[['Neighbor', 'Alt', "Seq"]]


def plot96Mtrx(df, height='Seq', neighbor='Neighbor', mut='Alt', title='', rows=None):
    ''' Trinucleotide Plot
  
    Parameters
    ----------
    df : pandas.DataFrame
       Frequency table, which can be obtained by genTrinucleotideMtrx.
        - columns names: [height, mut, neighbor] 
          - Neigbor: Nucleotide nearby snv. e.g A_C 
          - mut: SNV
          - height: frequency of mutation grouped by mut and neighbor
        - index: SNV ["T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G"]  

    height : str, optional
       Frequency of mutation grouped by mut and neighbor (the default is 'Seq')
    neighbor : str, optional
       Nucleotide nearby snv. e.g. A_C(the default is 'Neighbor')
    mut : str, optional
       SNV (the default is 'Alt')
    title : str, optional
       Figure title (the default is '')
    rows : str, optional
       Factor that group figures.(the default is None)
  
    Returns
    -------
    [type]
      [description]
    '''

    tri_n = OrderedDict([('A_A', 0), ('A_C', 0), ('A_G', 0), ('A_T', 0),
                         ('C_A', 0), ('C_C', 0), ('C_G', 0), ('C_T', 0),
                         ('G_A', 0), ('G_C', 0), ('G_G', 0), ('G_T', 0),
                         ('T_A', 0), ('T_C', 0), ('T_G', 0), ('T_T', 0),
    ])

    mut_panel = ['C>A', 'C>G', 'C>T', 'A>G', 'A>C', 'A>T']
    colors = ["deepskyblue", "black", "red",
              "lightgray", 'springgreen', "pink"]
    if not rows is None:
        row_panel = df[rows].unique()
    else:
        row_panel = [title]
    num_row = len(row_panel)
    fig, axarr = plt.subplots(num_row, 6, sharex=True,
                              sharey=True, figsize=(24, 2*num_row))
    #LEN: HACK- if there's just one plot, i.e. no TCGA panel,
    #then axarr = 2-D np-array with just one elm
    if num_row == 1:
        axarr = np.array([axarr])

    for j, r in enumerate(row_panel):
        for i, c in enumerate(mut_panel):
            tmp_tri = pd.Series(tri_n)
            determiner = (df[mut] == c)
            if num_row > 1:
                determiner = determiner & (df[rows] == r)

            neighbor_h = pd.Series(
                df.loc[determiner, [neighbor, height]].set_index(neighbor).to_dict()[height])
            neighbor_h /= neighbor_h.sum()
            tmp_tri.update(neighbor_h)
            ax = axarr[j, i]
            if j == 0:
                ax.set_title(c, weight='bold')
                ax.spines['top'].set_visible(False)
            if i < len(mut_panel) - 1:
                ax.spines['right'].set_visible(False)

            if i == 0 and num_row == 1:
                ax.set_ylabel('Relative Contribution', weight="bold")
            if i == 0 and num_row > 1:
                ax.set_ylabel(r, weight="bold", fontsize=12)
            tmp_tri.plot(kind='bar', ax=ax, color=colors[i], width=.9)
    fig.subplots_adjust(wspace=0, hspace=0.2)
    return fig


def convertTriMatrix(df, height='Seq', neighbor='Neighbor', mut='Alt'):
    '''Convert the Trtrinucleotide matrix from the pandas dataframe object 
    created by genTrinucleotideMtrx into the following dictionary:
    {‘CtoA’: {‘a_a’: 5, ‘a_c’:10, ‘a_g’: 0, …},
     ‘CtoG’: {‘a_a’: 5, ‘a_c’:10, ‘a_g’: 0, …},
    ...
    }
    Parameters
    ----------
    df : pandas.DataFrame
       Frequency table, which can be obtained by genTrinucleotideMtrx.
        - columns names: [height, mut, neighbor] 
          - Neigbor: Nucleotide nearby snv. e.g A_C 
          - mut: SNV
          - height: frequency of mutation grouped by mut and neighbor
        - index: SNV ["T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G"]  
    Returns
    -------
    dictionary where keys = ['C>A', 'C>G', 'C>T', 'A>G', 'A>C', 'A>T'] and the
    values of each key is a dictionary of neighbor counts, e.g.
    {‘a_a’: 5, ‘a_c’:10, ‘a_g’: 0, …}
    '''
    tri_n = OrderedDict([('A_A', 0), ('A_C', 0), ('A_G', 0), ('A_T', 0),
                         ('C_A', 0), ('C_C', 0), ('C_G', 0), ('C_T', 0),
                         ('G_A', 0), ('G_C', 0), ('G_G', 0), ('G_T', 0),
                         ('T_A', 0), ('T_C', 0), ('T_G', 0), ('T_T', 0),
    ])

    mut_panel = ['C>A', 'C>G', 'C>T', 'A>G', 'A>C', 'A>T']

    ret = {}
    for c in mut_panel:
        tmp_tri = pd.Series(tri_n)
        determiner = (df[mut] == c)
        #Get the un-normalized counts
        neighbor_h = pd.Series(
            df.loc[determiner, [neighbor, height]].set_index(neighbor).to_dict()[height])
        tmp_tri.update(neighbor_h)

        #Integrate them into the overall dictionary 'ret'
        ret[c] = dict([(i, int(tmp_tri[i])) for i in tmp_tri.index])
        #print(ret)
        #sys.exit()
    #HACK to match David's convention: instead of 'A>G', 'A>C', 'A>T'
    #David is using T>C, T>G, T>A
    conv = {'A>G':'T>C', 'A>C': 'T>G', 'A>T':'T>A'}
    for c in conv:
        new = conv[c]
        ret[new] = ret.pop(c)
    
    return ret

def main(sampleID, output,maf,ref_fasta,ref_cancer,json_output):
    tri_mtrx = genTrinucleotideMtrx(maf=maf, ref_fasta = ref_fasta)
    if ref_cancer and len(ref_cancer): #non-empty list of TCGA refs
        tri_mtrx['Group'] = sampleID
        #make list of TCGA ref pd matrices
        mtrxs = [pd.read_table(ref_f) for ref_f in ref_cancer]
        mtrxs.insert(0,tri_mtrx)
        fig = plot96Mtrx(df=pd.concat(mtrxs), rows='Group')
    else:
        fig = plot96Mtrx(df=tri_mtrx)
    fig.savefig('{}.pdf'.format(output),dpi=200)
    plt.clf()

    #Save tri_mtrx as json--need to convert dataFrame to this format:
    #{‘C>A’: {‘a_a’: 5, ‘a_c’:10, ‘a_g’: 0, …},
    # ‘C>G’: {‘a_a’: 5, ‘a_c’:10, ‘a_g’: 0, …},
    #}
    tri_mtrx_dict = convertTriMatrix(tri_mtrx)
    with open(json_output, 'w') as f:
        json.dump(tri_mtrx_dict, f, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Tools to download public genomic data", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', type=str,
                        required=True, help="Output Path")
    parser.add_argument('-j', '--json_output', type=str,
                        required=True, help="Json Output Path")
    parser.add_argument('-r', '--ref', type=str,
                        required=True, help="Reference Fasta file")
    parser.add_argument('-m', '--maf', type=str,
                        required=True, help="MAF File")

    parser.add_argument('-c', '--cancer', type=str, action='append',
                        required=False, help="Reference cancer mutation frequency file")
    parser.add_argument('-n', '--name', type=str, default='sample', required=False, help="Sample Name")
  

    args = parser.parse_args()
    main(maf=args.maf, ref_fasta=args.ref, output=args.output, ref_cancer=args.cancer, sampleID=args.name, json_output=args.json_output)
