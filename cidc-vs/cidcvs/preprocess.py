#!/usr/bin/env python3

import os,pybedtools
from io import StringIO
import pandas as pd
import numpy as np
from .definition import *

def addMeta(df,name_list):
    meta_info = [x.split('.')[0] for x in name_list]
    df['Center'] = [x.split('_')[0] for x in meta_info]
    df['Source'] = [x.split('_')[1] for x in meta_info]
    df['ID'] = [x.split('_')[2] for x in meta_info]

    return df

def concatPurity(base_path):
    file_list = os.listdir(base_path)
    df = pd.concat(
        [pd.read_table('{0}/{1}'.format(base_path, x)) for x in file_list],
         axis=0)

    # Add metadata
    df = addMeta(df=df,name_list=file_list)

    return df


def concatClonality(base_path):
    file_list = os.listdir(base_path)
    df = pd.DataFrame({'Clonality':[ estClonality( 
                                  pd.read_table('{0}/{1}'.format(base_path,x))
                                  ) for x in file_list],
                       })

    # Add metadata
    df = addMeta(df=df, name_list=file_list)

    return df

def concatCNV(base_path):
    file_list = os.listdir(base_path)
    df = pd.concat([pd.read_table('{0}/{1}'.format(base_path, x),
                                   index_col=0).drop(['Gene ID', 'Cytoband'], axis=1) for x in file_list],
                    axis=1,
                    join='inner'
                   )
    cor_df = df.corr()

    # format metadata
    metadata = pd.DataFrame({
                            'Center': df.columns.map(lambda x: x.split('_')[0]),
                            'Source': df.columns.map(lambda x: x.split('_')[1]),
                            'ID': df.columns.map(lambda x: x.split('_')[2])
                              },index=df.columns)

    return cor_df, metadata[['Center','Source','ID']]


def fetchMut(maf, only_protein=True, dna_alt_col="HGVSc", mut_type = 'snv'):
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
        df = df.loc[ (df["Variant_Type"] == "SNP") & ( ~df[dna_alt_col].isna()), :]
    elif mut_type == 'indel':
        df = df.loc[(df["Variant_Type"].isin(['INS','DEL'])) &
                    (~df[dna_alt_col].isna()), :]
    else:
        raise KeyError('{} is not valid, Only accept: indel, snv.'.format(mut_type))
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


def extractMutCol(maf, dna_alt_col="HGVSc", **kwargs):
    filter_maf = fetchMut(maf, dna_alt_col=dna_alt_col, **kwargs)
    return filter_maf["Hugo_Symbol"]+'_'+filter_maf[dna_alt_col]


def calJaccardMtrx(base_path, **kwargs):
    files_list = os.listdir(base_path)
    mut_list = [ (x.split('.')[0],
                  extractMutCol("{0}/{1}".format(base_path, x), **kwargs)
                  ) for x in files_list
                ]
    df = []
    for r in mut_list:
        for c in mut_list:
            df.append([r[0],c[0],jaccard(r[1],c[1])])

    df = pd.DataFrame(df,columns=['rows','cols','Jaccard Index'])
    df = pd.pivot_table(df, index='rows', columns='cols',
                        values='Jaccard Index')
    metadata = pd.Series(list(df.columns.map(lambda x: x.split('_')[0])),
                        index=df.columns)
 
    return df, metadata

def genTrinucleotideMtrx(maf, ref_fasta, dna_alt_col="HGVSc"):
    conv = dict(zip(('A>G', 'T>C', 'C>T', 'G>A', 'A>T', 'T>A', 'A>C', 'T>G', 'C>A', 'G>T', 'C>G', 'G>C'),
                    ("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G")))

    snv_maf = fetchMut(maf, dna_alt_col=dna_alt_col).iloc[:, 1:]

    snv_intervel = pybedtools.BedTool.from_dataframe(snv_maf)
    ref_fasta = pybedtools.example_filename(ref_fasta)
    snv_obj = snv_intervel.sequence(fi=ref_fasta, tab=True, name=True)
    snv_seq = StringIO(open(snv_obj.seqfn).read())

    df = pd.read_table(snv_seq, names=[dna_alt_col, 'Seq'])
    df['Neighbor'] = df.Seq.map(lambda x: x[0]+"_"+x[2])
    df['Alt'] = df[dna_alt_col].map(lambda x: x[-3:]).map(conv)
    return df.groupby(['Neighbor', 'Alt']).count().reset_index()[['Neighbor', 'Alt', "Seq"]]

def iterMaf(base_path,ref_fasta):
    files_list = os.listdir(base_path)
    for x in files_list:
        yield (x.split('.')[0],
               genTrinucleotideMtrx(
                   maf="{0}/{1}".format(base_path, x), ref_fasta=ref_fasta)
               )
   
