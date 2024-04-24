#single sync check
#expected output to STDOUT:
#Five numbers, corresponding to (and in this order):
# Number of minor alleles that are exactly one read depth
# Number of minor alleles that are >1 but fail minor threshold
# Number of sites with 3 alleles
# Number of 3 allele sites that are exaclty one read depth
# Average depth of all *considered* sites (not total avg depth)


import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Compare to sync files")
parser.add_argument('--p1', help='Parental Sync File One', required=True, type=str)
parser.add_argument('--rc', help='Minimum read count in Parental sync files (default=12)', required=False, default=12, type=int)
parser.add_argument('--mc', help='Minimum minor read count (default=1)', required=False, default=1, type=int)
parser.add_argument('--mt', help='Parental minor count minimum read count frequency to major (default=0.25[25%])', required=False, default=0.25, type=float)
#parser.add_argument('--o', help='Output File prefix (default=mismatch_check_results.csv)', required=False, default="mismatch_check_results.csv", type=str)
args = parser.parse_args()

### Load sync file and assign column names
df = pd.read_csv(args.p1, sep='\t', header=None) # parent 1
df.columns = ['chr', 'pos', 'refA', 'countA']
df[['A1', 'T1', 'C1', 'G1', 'countN1', 'countDel1']] = df['countA'].str.split(':', expand=True)

###Drop read columns with less than the minimum read count
df = df[(df['A1'].astype(int) + df['T1'].astype(int) + df['C1'].astype(int) + df['G1'].astype(int)).ge(args.rc)]

allele_dict = {'0':'A','1':'T','2':'C','3':'G'}
### Get Major/Minor
def major_minor(row, set_num):
	if set_num == 1:
		A_count = int(row['A1'])
		T_count = int(row['T1'])
		C_count = int(row['C1'])
		G_count = int(row['G1'])
	else:
		sys.exit("You just got wrecked")
	
	counts = np.array([A_count, T_count, C_count, G_count])
	count_index = np.argsort(counts)
	major_allele = allele_dict.get(str(count_index[3])) #Major count (Base with greatest read count)
	
	if counts[count_index[2]] == 0:
		minor_allele = "N"
	else:
		minor_allele = allele_dict.get(str(count_index[2]))
	if counts[count_index[1]] == 0:
		third_allele = "N"
	else:
		third_allele = allele_dict.get(str(count_index[1]))
	if counts[count_index[0]] == 0:
		fourth_allele = "N"
	else:
		fourth_allele = allele_dict.get(str(count_index[0]))

	if counts[count_index[2]] == 1:
		singleton = 1
	else:
		singleton = 0

	if (counts[int(count_index[2])] > args.mc) and (counts[int(count_index[2])] <= (float(counts[int(count_index[3])])* args.mt)):
		minorfail = 1
	else:
		minorfail = 0

	if counts[count_index[1]] >= 1:
		third = 1
	else:
		third = 0

	if counts[count_index[1]] == 1:
		thirdcount = 1
	else:
		 thirdcount = 0
	
	depth = int(A_count) + int(T_count) + int(G_count) + int(C_count)
	
	return singleton, minorfail, third, thirdcount, depth

### Run major_minor function to get sync1 and then sync2 counts
try:
	df[["single1", "s1", "s2", "s3", "d"]] = df.apply(lambda row: major_minor(row,1), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	sys.exit("Empty Dataframe (Error catch 1) Exiting....")

total1=df['single1'].astype(int).sum()
total2=df['s1'].astype(int).sum()
total3=df['s2'].astype(int).sum()
total4=df['s3'].astype(int).sum()
total5=(df['d'].astype(int).sum()/(df.shape[0]-1))
print(str(total1))
print(str(total2))
print(str(total3))
print(str(total4))
print(str(total5))

#optional print final df that lists M and MM columns by site (sanity check purposes)
#df.to_csv('test.csv', sep=',', index=False, header=True)
