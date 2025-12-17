# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 08:35:34 2022

@author: 123
"""

import pandas as pd
from tidepy.pred import TIDE
import os
os.chdir('E:/work/240624_brca/27_tide')

df=pd.read_csv("data/2clusters_gsva_nes.tsv",sep="\t")
import math
df1=df-df.mean()
result=TIDE(df1, cancer='Other',pretreat=False,vthres=0.)
result.to_csv('result/brca_train_tide.tsv',sep="\t")


# [WARN] 4.46 % Genes are missing after converting to Entrez ID
# [WARN] 13.4% MSI signature genes are missing on input expression profile.