'''
Summarize the expression level from the output (transcript.featureCounts) of Tao's RNAseq pipeline, output the table with one column per sample and one row per gene. 

The example to run this script:

- for human samples
python3 summarizeByCoverage.py --designFile human.design4expressionsum.txt --countsFileName transcript.featureCounts --refFlat hs38d1.fa_cnvkit.txt --sumFile human.sum.txt 

- for mouse samples
python3 summarizeByCoverage.py --designFile mouse.design4expression.txt.sum --countsFileName transcript.featureCounts --refFlat mm10.fasta_cnvkit.txt --sumFile human.sum.txt 

To get help for each parameter of the script, please run: python3 summarizeByCoverage.py -h
'''


import time
import os
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


# In refFlat file, one gene symbol (name) may correspond to multiple transcript IDs (gene id from transcript.featureCounts file)
# geneid2gene: {geneid:gene, ...}
def getGeneid2gene(refFlat):
	csvfile = refFlat
	geneid2gene = {}
	rows = tools.readCsvFile(csvfile, delimiter='\t')
	for row in rows:
		gene = row[0].strip()
		geneid = row[1].strip().split('.', maxsplit=1)[0]
		geneid2gene[geneid] = gene
	return geneid2gene

# average the expression level of mutliple transcripts per a gene, namely, expression level of a gene
def avgTranscriptsPerGene(geneDic):
	geneDicAvg = {}
	for gene, samples in geneDic.items():
		geneDicAvg[gene] = {sample:0 for sample in samples}
		for sample in samples:
			transcripts = geneDic[gene][sample]
			geneDicAvg[gene][sample] = sum([transcript for transcript in transcripts]) / len(transcripts)
	return geneDicAvg

# Return geneDicnew
# geneDicnew: {gene:samples, ...}
# sample: {sample:transcripts, ...}
# transcripts: [exp, ...], exp is the expression level (counts) of the specific transcript
def replaceGeneidByGene(geneDic, geneid2gene):
	geneDicnew = {}

	# addup in each gene the expression level values of all transcripts from the same gene
	for geneid,samples in geneDic.items():
		if geneid in geneid2gene.keys():
			gene = geneid2gene[geneid]
			if gene not in geneDicnew.keys():
				geneDicnew[gene] = {}
			for sample in samples:
				if sample not in geneDicnew[gene].keys():
					geneDicnew[gene][sample] = []
				geneDicnew[gene][sample].append(geneDic[geneid][sample])
		else:
			print('Warning, transcript ID (Geneid) {} from transcript.featureCounts files is not in refFlat file'.format(geneid))
			# keep the geneid in geneDic if geneid is not in geneid2gene (refFlat file)
			#geneDicnew[geneid] = geneDic[geneid]
	
	return geneDicnew

# Return expressions
# expressions: {sample:rows, ...}
# rows: [dic, ...]
# dic: {'Geneid':Geneid, ..., 'expression':expression}
def getExpressions(margs, samples):
	expressions = {}
	for args,sample in zip(margs, samples):
		tumor, path2output = args
		csvfile = path2output
		rows = tools.readDictCsvFile(csvfile, delimiter='\t')
		rowsNew = []
		for row in rows:
			rowNew = {'Geneid':row['Geneid'].strip(), 
				'expression':float(row['expression'].strip())}
			rowsNew.append(rowNew)
		expressions[sample] = rowsNew
	return expressions

# return geneDic
# geneDic: {geneid: counts4samples, ...}
# counts4samples: {sample: counts, ...}
#
# expressions: {sample:counts, ...}
# counts: [..., rpkm]
def summarize4expressionGeneLevel(expressions):
	geneDic = {}
	for sample,rows in expressions.items():
		for row in rows:
			gene = row['Geneid'].split('.', maxsplit=1)[0]
			counts = row['expression']
			if gene not in geneDic.keys():
				geneDic[gene] = {}
			if sample not in geneDic[gene].keys():
				geneDic[gene][sample] = 0
			geneDic[gene][sample] += counts

	return geneDic


# add qc info from csvfile to qcDic
def getqc(qcDic, sampleid, csvfile, seqid=''):
	rows4qc = tools.readCsvFile(csvfile, delimiter='\t')
	for row4qc in rows4qc:
		qcResult = row4qc[0].strip()
		qc = row4qc[1].strip()
		seqfile = row4qc[2].strip()
		qcid = (qc, seqid)
		qcid = '/'.join(qcid)
		if qcid not in qcDic.keys():
			qcDic[qcid] = {}
		qcDic[qcid][sampleid] = qcResult

# qcDic: {qcid:sample, ...}
# qcid: character string, 'qc, seqfile'
# qc: character string, e.g. 'Basic Statistics', 'Per base sequence quality', 'Per tile sequence quality', 'Per sequence quality scores'
# seqfile: character string, e.g. 'ZHE718-27-5985_T_RNA_wholernaseq-24-13601_S1_R1_001.fastq.gz'
# sample: {sampleid:qcResult, ...}
# sampleid: character string, e.g. 'IL2_104M' 
# qcResult: character string, e.g. 'PASS', 'FAIL', 'WARN'
def summarize4qc(designRows):
	qc1 = "qc_summary_1.txt"
	qc2 = "qc_summary_2.txt"
	qcDic = {}
	sampleids = []
	for row in designRows:
		sampleid = row['Tumor_sample_ID'].strip()
		csvfile = os.path.join(row['Path2output'].strip(), qc1)
		if os.path.isfile(csvfile):
			getqc(qcDic, sampleid, csvfile, seqid=qc1)
		else:
			print('{} not found, ignore it in summarization of qc_summary'.format(csvfile))
		csvfile = os.path.join(row['Path2output'].strip(), qc2)
		if os.path.isfile(csvfile):
			getqc(qcDic, sampleid, csvfile, seqid=qc2)
		else:
			print('{} not found, ignore it in summarization of qc_summary'.format(csvfile))
		sampleids.append(sampleid)
	return (qcDic, sampleids)


def mainFunc(args4main):
	print('start', datetime.datetime.now())
	designFile = args4main['designFile']
	designRows = tools.readDictCsvFile(designFile, delimiter='\t')

	countsFileName = args4main['countsFileName']
	margs = []
	samples = []
	for row in designRows:
		#normal = row['Normal_sample_ID'].strip()
		tumor = row['Tumor_sample_ID'].strip()
		path2output = os.path.join(row['Path2output'], countsFileName)
		args = (tumor, path2output)

		margs.append(args)
		samples.append(tumor)

	print("Start reading featureCounts file at", datetime.datetime.now())
	# get expression level values for the specific transcript (Geneid in featureCounts file) in the specific sample
	expressions = getExpressions(margs, samples)

	geneDic = summarize4expressionGeneLevel(expressions)

	# Map all transcripts of the same gene to the gene name
	refFlat = args4main['refFlat']
	print("Start reading refFlat ({}) at {}".format(refFlat, datetime.datetime.now()))
	geneid2gene = getGeneid2gene(refFlat)

	# addup the expression level values of all transcripts from the specific gene
	print("Start replacing gene ID with gene name at", datetime.datetime.now())
	#print('before replace geneid, geneDic[NR_075077][Double15T1]', geneDic['NR_075077']['Double15T1'])
	geneDic = replaceGeneidByGene(geneDic, geneid2gene)
	#print('after replace geneid, geneDic[C1orf141][Double15T1]', geneDic['C1orf141']['Double15T1'])

	# average the expression level of mutliple transcripts per a gene, namely, expression level of a gene
	geneDicAvg = avgTranscriptsPerGene(geneDic)
	#print('after average, geneDic[C1orf141][Double15T1]', geneDicAvg['C1orf141']['Double15T1'])

	sumFile = args4main['sumFile'].strip()
	geneDic = geneDicAvg
	tools.outputSumFile(geneDic, sumFile, samples)

	# summarize the qc of raw sequences
	qcDic, samples = summarize4qc(designRows)
	tools.addSamples(qcDic, samples)
	sumFile4qc = sumFile + '.qc'
	tools.outputSumFile(qcDic, sumFile4qc, samples)
	#
	print('end', datetime.datetime.now())


if __name__ == "__main__":
	import argparse
	import textwrap

	# Parse command line arguments
	descriptStr = '''Process all samples in a design file and produce the summarization about the gene expression for samples'''
	parser = argparse.ArgumentParser(prog='summarizeByCoverage.py', description = textwrap.dedent(descriptStr), 
			formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument(
		'--designFile', 
		required = True,
		type=str,
		default='',
		help = "Design file with header line. It contains three columns separated by tab: Patient_ID, \
			Tumor_sample_ID, Path2output. Patient_ID is identifier of patient, Tumor_sample_ID \
			is the identifer of tumor sample, Path2output is the output directory of Tao's RNAseq pipeline \
			for the tumor sample. The column name in header line must be Patient_ID, \
			Tumor_sample_ID, Path2output.",
		)

	parser.add_argument(
		'--countsFileName', 
		required = True,
		type=str,
		default='transcript.featureCounts',
		help = 'Feature count file of transcripts, e.g. transcript.featureCounts',
		)

	parser.add_argument(
		'--refFlat', 
		required = True,
		type=str,
		default='',
		help = 'RefFlat file mapping gene symbol to transcript ID, e.g. \
			/project/shared/xiao_wang/data/hg38/hs38d1.fa_cnvkit.txt for human samples, \
		 	/project/shared/xiao_wang/data/mm10/mm10.fasta_cnvkit.txt for mouse samples',
		)

	parser.add_argument(
		'--sumFile', 
		required = True,
		type=str,
		default='',
		help = 'Output file, which is a csv like file with columns separated by tab',
		)

	args = parser.parse_args()

	args4mainFunc =	{
				'designFile': args.designFile.strip(),
				'countsFileName': args.countsFileName.strip(),
				'refFlat': args.refFlat.strip(),
				'sumFile': args.sumFile.strip(),
			}

	mainFunc(args4mainFunc)

