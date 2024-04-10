#More ancestry information for F2 
#expected output, one text file with a single chrom,match,mismatch row
#MJL 04/02/24

import argparse
import numpy as np
import pandas as pd
import sys
import os.path
import csv

parser = argparse.ArgumentParser(description="Compare to sync files")
parser.add_argument('--p1', help='Parental Sync File One', required=True, type=str)
parser.add_argument('--p2', help='Parental Sync File Two', required=True, type=str)
parser.add_argument('--f', help='Input F2 Sync File', required=True, type=str)
parser.add_argument('--rc', help='Minimum read count in Parental sync files (default=12)', required=False, default=12, type=int)
parser.add_argument('--mc', help='Minimum minor read count (default=1)', required=False, default=1, type=int)
parser.add_argument('--mt', help='Parental minor count minimum read count frequency to major (default=0.25[25%])', required=False, default=0.25, type=float)
#parser.add_argument('--o', help='Output File prefix (default=mismatch_check_results.csv)', required=False, default="mismatch_check_results.csv", type=str)
args = parser.parse_args()

### Load sync files and assign column names
df_s1 = pd.read_csv(args.p1, sep='\t', header=None) # parent 1
df_s2 = pd.read_csv(args.p2, sep='\t', header=None) # parent 2
df_s3 = pd.read_csv(args.f, sep='\t', header=None) # parent 3
df_s1.columns = ['chr', 'pos', 'refA', 'countA']
df_s2.columns = ['chr', 'pos', 'refA', 'countB']
df_s3.columns = ['chr', 'pos', 'refA', 'countC']

### Merge Sync files at common sites
df = pd.merge(df_s1[['pos', 'countA']], df_s2[['pos', 'countB']], on='pos')
df[['A1', 'T1', 'C1', 'G1', 'countN1', 'countDel1']] = df['countA'].str.split(':', expand=True)
df[['A2', 'T2', 'C2', 'G2', 'countN2', 'countDel2']] = df['countB'].str.split(':', expand=True)
df.drop(df.columns[[1,2,7,8,13,14]], axis=1, inplace=True) #drop unneeded columns(total counts, Indel and Del counts)

###Drop read columns with less than the minimum read count
df = df[(df['A1'].astype(int) + df['T1'].astype(int) + df['C1'].astype(int) + df['G1'].astype(int)).ge(args.rc)]
df = df[(df['A2'].astype(int) + df['T2'].astype(int) + df['C2'].astype(int) + df['G2'].astype(int)).ge(args.rc)]

allele_dict = {'0':'A','1':'T','2':'C','3':'G'}
### Get Major/Minor
def major_minor(row, set_num):
	if set_num == 1:
		A_count = int(row['A1'])
		T_count = int(row['T1'])
		C_count = int(row['C1'])
		G_count = int(row['G1'])
	elif set_num == 2:
		A_count = int(row['A2'])
		T_count = int(row['T2'])
		C_count = int(row['C2'])
		G_count = int(row['G2'])
	else:
		print("Sanity Error Print (function majorminor_grouped)")
	counts = np.array([A_count, T_count, C_count, G_count])
	count_index = np.argsort(counts)
	major_allele = allele_dict.get(str(count_index[3])) #Major count (Base with greatest read count)
	### Minor allele count threshold and percent of major threshold
	if (counts[int(count_index[2])] >= args.mc) and (counts[int(count_index[2])] >= (float(counts[int(count_index[3])])* args.mt)):
		minor_allele = allele_dict.get(str(count_index[2]))
	else:
		minor_allele = 'N'

	return major_allele, minor_allele

### Run major_minor function to get sync1 and then sync2 counts
try:
	df[["1A1", "1A2"]] = df.apply(lambda row: major_minor(row, 1), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	sys.exit("Empty Dataframe (Error catch 1) Exiting....")
try:
	df[["2A1", "2A2"]] = df.apply(lambda row: major_minor(row, 2), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	sys.exit("Empty Dataframe (Error catch 1) Exiting....")


#Merge primary dataframe to offspring sync file on site, drop unneeded cols, split count column
df.drop(df.columns[[1,2,3,4,5,6,7,8]], axis=1, inplace=True) #drop unneeded columns
df_2 = pd.merge(df[['pos', '1A1', '1A2', '2A1', '2A2']], df_s3[['pos', 'countC']], on='pos')
df_2[['A', 'T', 'C', 'G', 'countN1', 'countDel1']] = df_2['countC'].str.split(':', expand=True)
df_2.drop(df_2.columns[[5,10,11]], axis=1, inplace=True) #drop unneeded columns(total counts, Indel and Del counts)

def mismatcher(row):
	match_count = 0
	mismatch_count = 0
	#Get bases represented in parental syncs
	group1 = [str(row['1A1']), str(row['1A2'])]
	group2 = [str(row['2A1']), str(row['2A2'])]
	alleles = list(np.unique(group1 + group2))
	#Remove 'N' values if present
	try:
		alleles.remove('N')
	except ValueError:
		pass

	#Count matches
	for base in alleles:
		match_count += int(row[base])
	#Count mismatch (all reads not in overlap)
	total_reads = int(row['A']) + int(row['T']) + int(row['C']) + int(row['G'])
	mismatch_count = total_reads - match_count

	return match_count, mismatch_count

try:
	df_2[["M", "MM"]] = df_2.apply(lambda row: mismatcher(row), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	sys.exit("Empty Dataframe (Error catch 2) Exiting....")

match_total = df_2['M'].sum()
mismatch_total = df_2['MM'].sum()

chrom_dict = {'3':'2L','4':'X','5':'3L','7':'2R', '8':'3R'}
out_ind = args.f.split('.')[1]
chrom_index = args.f.split('.')[0]
out_chrom = chrom_dict.get(str(chrom_index))

outname = "mismatch_check_results.csv" #Match initited file
exist_check = os.path.isfile(outname)
with open (outname, 'a') as outfile:
	topline = ['ind','chrom','match','mismatch']
	writer = csv.DictWriter(outfile, delimiter=',', lineterminator='\n',fieldnames=topline)
	if not exist_check:
		writer.writeheader()
	writer.writerow({'ind':str(out_ind),'chrom':str(out_chrom),'match':str(match_total), 'mismatch':str(mismatch_total)})

#optional print final df that lists M and MM columns by site (sanity check purposes)
#df_2.to_csv('test.csv', sep=',', index=False, header=False)
