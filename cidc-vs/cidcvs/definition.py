#!/usr/bin/env python3
import numpy as np

def estClonality(df):
    P = df['size'] / df['size'].sum()
    n = df.shape[0]
    return 1 + (P*np.log2(P)).sum()/np.log2(n)


def jaccard(a, b):
    a, b = set(a), set(b)
    union = a.intersection(b)
    return float(len(union)) / (len(a) + len(b) - len(union))


        
   
