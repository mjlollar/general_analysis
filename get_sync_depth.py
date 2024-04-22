#Average Depth of sync file
#MJL 04/22/24

import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Compare to sync files")
parser.add_argument('--s', help='Sync File', required=True, type=str)
args = parser.parse_args()

### Load sync file and assign column names
df = pd.read_csv(args.s, sep='\t', header=None)
df.columns = ['chr', 'pos', 'refA', 'countA']
df[['A', 'T', 'C', 'G', 'countN1', 'countDel1']] = df['countA'].str.split(':', expand=True)

read_list = ['A', 'T', 'C', 'G']
df['depth'] = df[read_list].astype(int).sum(axis=1)
final_depth = (df['depth'].sum()/(df.shape[0]-1))
print(str(final_depth))
