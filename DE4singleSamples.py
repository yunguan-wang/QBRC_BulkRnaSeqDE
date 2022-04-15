'''
Calculate differential expression levels between genes from two samples using the equations below:

a=log2(sample1_counts_each_gene)-log2(average(sample1_counts_all_genes))

b=log2(sample2_counts_each_gene)-log2(average(sample2_counts_all_genes))

log2FoldChange = b - a


Input:
gene-sample raw counts: gene-sample table produced by summarize4expression.py, the gene expression level in table with one column per sample and one row per gene. 

group design file: same as the group design file required by DE.r

Output:
differential expression of genes of control vs sample: 2 columns, column 1 is gene symbol, column 2 is fold change.


The example to run this script:

- for human samples
python3 DE4singleSamples.py --group groups.txt --geneSampleCounts 155.193.237.v3.counts.rev_stranded.csv  --sumFile DE4singleSamples.csv


To get help for each parameter of the script, please run: python3 DE4singleSamples.py -h
'''


import time
import os
import math
import datetime


# Import tools.py in different cases
import pathlib

rootDir = os.path.dirname(__file__)
fullPath2tools = os.path.join(rootDir, "tools.py")
if os.path.isfile(fullPath2tools):
	import tools
elif os.path.isfile(os.path.join(rootDir, "script", "tools.py")):
	# assume tools.py in script subdirectory
	from script import tools


# convert dicList to dicDic
# dicList: [dic, ...]
# dic: {k:v, ...}
#
# dicDic: [key:dic, ...]
# dic: {k:v, ...}
def dicList2dicDic(geneSamples):
	dicDic = {}
	for row in geneSamples:
		gene = row['']
		for sample, v in row.items():
			if sample not in dicDic.keys():
				dicDic[gene] = {}
			dicDic[gene][sample] = v

# Return normalized log2 of counts for each gene in each sample
#
# Input
# geneSamples: [samples4gene, ...]
# samples4gene: {'':genename, sample:counts, ...}
#
def normByAvg(geneSamples):
	ngenes = len(geneSamples)
	sum4sample = {k:0 for k,v in geneSamples[0].items() if k != ''}

	for row in geneSamples:
		for sample,counts in row.items():
		# dic: {'':gene, sample:counts, ...}
			if sample == '':
				continue
			sum4sample[sample] += counts
	log2avg = {}
	for sample,c in sum4sample.items():
		avg = c/ngenes
		# set avg=1 if avg is a flost number close to zero 
		if not (avg > 0):
			avg = 1
		log2avg[sample] = math.log2(avg)
			
	# normalization
	norms4geneSamples = {}
	for row in geneSamples:
		gene = row['']
		norm = {}
		for k,v in row.items():
			if k == '':
				continue
			if not (v > 0):
				v = 1
			norm[k] = math.log2(v) - log2avg[k]
		norms4geneSamples[gene] = norm
	return norms4geneSamples
			
def str2int(geneSamples):
	diclist = []
	for row in geneSamples:
		dic = {}
		for k,v in row.items():
			if k == '':
				dic[k] = v
			else:
				dic[k] = int(v)
		diclist.append(dic)
	return diclist

# return (idpairs, id2group)
# idpairs: [pair, ...]
# pair: [ref, contrast]
def getpairs(rowsDesign):
	id2group = {}
	for row in rowsDesign:
		id2group[row['ID'].strip()] = row['Group'].strip()

	pairs = []
	for row in rowsDesign:
		referenceID = row['ID'].strip()
		contrasts = [item.strip() for item in row['Contrasts'].split(',')]
		for contrast in contrasts:
			pair = [referenceID, contrast]
			pairs.append(pair)

	idpairs = []
	for pair in pairs:
		for id,group in id2group.items():
			if pair[1] != group:
				continue
			idpairs.append([pair[0], id])
	return (idpairs, id2group)

# return foldchanges
# foldchanges: {gene:pair, ...}
# pair: {(refSample, contrastSample): foldchange, ...}
def compare(norms4geneSamples, pairs):
	foldchanges = {}
	for gene, samples in norms4geneSamples.items():
		foldchanges[gene] = {}
		for pair in pairs:
			pairkey = '_vs_'.join([pair[1], pair[0]])
			foldchanges[gene][pairkey] = samples[pair[1]] - samples[pair[0]]
	return foldchanges

# return dicDic but keys are replaced by the mapped group IDs
#  
def roottree2group(dicDic, id2group):
	dicDicNew = {}
	for gene,dic in dicDic.items():
		dicDicNew[gene] = {}
		for sample,value in dic.items():
			if sample in id2group.keys():
				sample = id2group[sample]
			elif '_vs_' in sample:
				s1,s2 = sample.split('_vs_')
				if s1 in id2group.keys():
					s1 = id2group[s1]
				if s2 in id2group.keys():
					s2 = id2group[s2]
				sample = '_vs_'.join([s2,s1])
			dicDicNew[gene][sample] = value
	return dicDicNew

def mainFunc(args4main):
	print('start at', datetime.datetime.now())

	group = args4main['group']
	rowsDesign = tools.readDictCsvFile(group, delimiter='\t')
	pairs, id2group = getpairs(rowsDesign)

	geneSampleCountsFile = args4main['geneSampleCounts']
	geneSamples = tools.readDictCsvFile(geneSampleCountsFile, delimiter='\t')
	geneSamples = str2int(geneSamples)


	# normalize the gene counts for each sample
	norms4geneSamples = normByAvg(geneSamples)

	# output normalized gene counts for each sample
	sumFile = args4main['sumFile'].strip()
	#norms4geneSamplesNew = roottree2group(norms4geneSamples, id2group)
	norms4geneSamplesNew = norms4geneSamples
	# get a gene name, any gene in the data table is OK.
	gene = next(iter(norms4geneSamplesNew))
	samplesNew = norms4geneSamplesNew[gene].keys()
	tools.outputSumFile(norms4geneSamplesNew, sumFile+'.csv', samplesNew, sep=',')

	# fold change between two samples
	foldchanges = compare(norms4geneSamples, pairs)

	# output fold changes for each sample pair
	foldchangeFile = sumFile + '.de.csv'
	#foldchangesNew = roottree2group(foldchanges, id2group)
	foldchangesNew = foldchanges
	# get a gene name, any gene in the data table is OK.
	gene = next(iter(foldchangesNew))
	pairkeys = foldchangesNew[gene].keys()
	tools.outputSumFile(foldchangesNew, foldchangeFile, pairkeys, sep=',')
	

	print('end at', datetime.datetime.now())


if __name__ == "__main__":
	import argparse
	import textwrap

	# Parse command line arguments
	descriptStr = '''Normalize the gene expression for each sample'''
	parser = argparse.ArgumentParser(prog='DE4singleSamples.py', description = textwrap.dedent(descriptStr), 
			formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument(
		'--group', 
		required = True,
		type=str,
		default='',
		help = "Design file with header line. Refer DE.r for details.",
		)

	parser.add_argument(
		'--geneSampleCounts', 
		required = True,
		default='',
		help = 'The gene-sample counts file produced by summarize4expression.py',
		)

	parser.add_argument(
		'--sumFile', 
		required = True,
		default='sum.csv',
		help = 'Output file which is a csv file with header, sum.csv by default',
		)

	args = parser.parse_args()

	args4mainFunc =	{
				'group': args.group.strip(),
				'geneSampleCounts': args.geneSampleCounts.strip(),
				'sumFile': args.sumFile.strip(),
			}

	mainFunc(args4mainFunc)

