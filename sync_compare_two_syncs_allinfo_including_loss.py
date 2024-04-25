#Comparing two sync files for %match/mismatch

## This version of sync comparison outputs additional detailed descriptions of 
## Major/Minor dynamics between sync files,
## like loss of minor allele in syncs, etc...

#MJL 03/20/24

import argparse
import datetime
import numpy as np
import pandas as pd
import sys

today = datetime.date.today()
parser = argparse.ArgumentParser(description="Compare to sync files")
parser.add_argument('--a', help='Sync File One', required=True, type=str)
parser.add_argument('--b', help='Sync File Two', required=True, type=str)
parser.add_argument('--rc', help='Minimum read count', required=True, type=int)
parser.add_argument('--mc', help='Minimum minor read count (default=2)', required=False, default=2, type=int)
parser.add_argument('--mt', help='Minor count minimum frequency (default=0.25[25%])', required=False, default=0.25, type=float)
parser.add_argument('--o', help='Output File prefix (default=date)', required=False, default=today.strftime("%b-%d-%Y"), type=str)
args = parser.parse_args()

### Load sync files and assign column names
df_s1 = pd.read_csv(args.a, sep='\t', header=None)
df_s2 = pd.read_csv(args.b, sep='\t', header=None)
df_s1.columns = ['chr', 'pos', 'refA', 'countA']
df_s2.columns = ['chr', 'pos', 'refA', 'countB']

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
	print("Getting Major/minor alleles....")
	#print("--------------------------------------------------------")
	#print("Listing any possible third allele detections (position, sync number):")
	df[["1A1", "1A2"]] = df.apply(lambda row: major_minor(row, 1), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	sys.exit("Empty Dataframe (Error catch 1) Exiting....")
try:
	df[["2A1", "2A2"]] = df.apply(lambda row: major_minor(row, 2), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	sys.exit("Empty Dataframe (Error catch 2), Exiting....")

#print("-------------------------------------------------------------")

### Compare Alleles
def allele_match(row):
	A1 = str(row['1A1'])
	a1 = str(row['1A2'])
	A2 = str(row['2A1'])
	a2 = str(row['2A2'])
	status1 = 'n'
	status2 = 'n'
	
	if A1 == A2: # Major alleles match 
		status1 = 'MAMA'
		if a1 == a2: # Minor alleles match (or N match)
			if a1 == 'N': # No minor in a1 and a2
				status2 = 'N'
			elif a1 != 'N': # True minor match
				status2 = 'M'
			else:
				sys.exit("Loop error code 1, exiting...")
		elif a1 != a2: # Minor mismatch (N in one line or 4 allelic site(true minor mismatch))
			if a1 == 'N': # Loss of minor allele in sync 1
				status2 = 'A'
			elif a1 != 'N':
				if a2 == 'N': # Loss of minor allele in sync 2
					status2 = 'B'
				elif a2 != 'N': # Both minor alleles present but do not match
					status2 = 'Q'
				else:
					sys.exit("Loop error code 2, exiting...")
			else: 
				sys.exit("Loop error code 3, exiting...")
		else:
			sys.exit("Loop error code 4, exiting...")
	
	elif A1 != A2: # Major alleles mismatch
		set_1 = [A1, a1]
		set_2 = [A2, a2]
		#remove N if present
		try:
			set_1.remove('N')
		except ValueError:
			pass
		try:
			set_2.remove('N')
		except ValueError:
			pass

		if not set(set_1).isdisjoint(set_2): #if there is an allele overlap
			status1 = "MAMI"
			if (len(set_1) == 2) and (len(set_2) == 2): #both minor alleles present
				if a1 == a2: #share minor alleles				
					status2 = "M"
				elif a1 != a2: #do not share minor alleles
					status2 = "Q"
				else:
					sys.exit("Loop error code 5, exiting...")
			elif (len(set_1) == 2) and (len(set_2) == 1): #minor present in sync1 only
				status2 = "B"
			elif (len(set_1) == 1) and (len(set_2) == 2): #minor present in sync2 only
				status2 = "A"
			else:
				sys.exit("Loop error code 6, exiting...")
		elif set(set_1).isdisjoint(set_2): #if there is no overlap of alleles
			status1 = "mismatch"
			if (len(set_1) == 2) and (len(set_2) == 2): #no overlap, two alleles each
				status2 = "QQ"
			elif (len(set_1) == 2) or (len(set_2) == 2): #minor present in one but no overlap
				status2 = "QQQ"
			elif (len(set_1) == 1) and (len(set_2) == 1): #only majors, no overlap
				status2 = "QQQQ"
			elif (len(set_1) == 1) or (len(set_1) == 1): #one major one minor, no overlap
				status2 = "QQQQQ"
			else:
				sys.exit("Loop error code 7, exiting...")
		else:
			sys.exit("Loop error code 8, exiting...")  
	else:
		sys.exit("Loop error code 9, exiting...")
	return status1, status2 

### Run allele matching function
try:
	print("Comparing alleles...")
	df[["status1", "status2"]] = df.apply(lambda row: allele_match(row), axis='columns', result_type='expand')
except ValueError: #catch instances where empty dataframe
	sys.exit("Empty Dataframe (Error catch 3), Exiting....")

### Tally status1
try: #Major/Major
	MAMA_tally = int(df['status1'].value_counts()['MAMA'])
except KeyError: #catch no counts of this value
	MAMA_tally = 0
try: #Major/Minor
	MAMI_tally = int(df['status1'].value_counts()['MAMI'])
except KeyError: #catch no counts of this value
	MAMI_tally = 0
try: #Mismatch
	mismatch_tally = int(df['status1'].value_counts()['mismatch'])
except KeyError: #catch no counts of this value
	mismatch_matches = 0

### Tally status2
try:
	N_tally = int(df['status2'].value_counts()['N'])
except KeyError:
	N_tally = 0
try:
	M_tally = int(df['status2'].value_counts()['M'])
except KeyError:
	M_tally = 0
try:
	A_tally = int(df['status2'].value_counts()['A'])
except KeyError:
	A_tally = 0
try:
	B_tally = int(df['status2'].value_counts()['B'])
except KeyError:
	B_tally = 0
try:
	Q_tally = int(df['status2'].value_counts()['Q'])
except KeyError:
	Q_tally = 0
try:
	QQ_tally = int(df['status2'].value_counts()['QQ'])
except KeyError:
	QQ_tally = 0
try:
	QQQ_tally = int(df['status2'].value_counts()['QQQ'])
except KeyError:
	QQQ_tally = 0
try:
	QQQQ_tally = int(df['status2'].value_counts()['QQQQ'])
except KeyError:
	QQQQ_tally = 0
try:
	QQQQQ_tally = int(df['status2'].value_counts()['QQQQQ'])
except KeyError:
	QQQQQ_tally = 0

#Some calcs
total_sites = MAMA_tally + MAMI_tally + mismatch_tally
overlap_tally = MAMA_tally + MAMI_tally
MAMA_per = (MAMA_tally / total_sites) * 100
MAMI_per = (MAMI_tally / total_sites) * 100
mismatch_per = (mismatch_tally / total_sites) * 100
overlap_per = (overlap_tally / total_sites) * 100
noloss_tally = N_tally + M_tally
total_sites_2 = A_tally + B_tally + noloss_tally
A_per = (A_tally / total_sites_2) * 100
B_per = (B_tally / total_sites_2) * 100
noloss_per = (noloss_tally / total_sites_2) * 100

### Write to output
outname = args.o + "_comparesyncresults.csv"
print("Writing to output, file name: " + outname)
outfile = open(outname, 'a')
outfile.write("Comparison between: " + str(args.a) + " and " + str(args.b) + "\n" )
outfile.write("description,raw_count,percentage_total" + "\n")
outfile.write("mismatch," + str(mismatch_tally) + "," + str(mismatch_per) + "\n")
outfile.write("overlap," + str(overlap_tally) + "," + str(overlap_per) + "\n")
outfile.write("Sync1_loss," + str(A_tally) + "," + str(mismatch_per) + "\n")
outfile.write("Sync2_loss," + str(B_tally) + "," + str(B_per) + "\n")
outfile.write("no_loss," + str(noloss_tally) + "," + str(noloss_per) + "\n")
outfile.write("total_sites_status1," + str(total_sites) + "\n")
outfile.write("total_sites_status2," + str(total_sites_2) + "\n")
outfile.write("Q_tally," + str(Q_tally) + "\n")
outfile.write("QQ_tally," + str(QQ_tally) + "\n")
outfile.write("QQQ_tally," + str(QQQ_tally) + "\n")
outfile.write("QQQ_tally," + str(QQQ_tally) + "\n")
outfile.write("QQQQ_tally," + str(QQQQ_tally) + "\n")
outfile.write("QQQQQ_tally," + str(QQQQQ_tally) + "\n")
outfile.close()


