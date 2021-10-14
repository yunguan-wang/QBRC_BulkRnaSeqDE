'''
Create two design files, one for job_expression.pl, another one for summarize4expression.py.

The example to run this script:

python3 design4expression.py example_data/manifest.txt result_path design_file.txt

Upon the command completing, two new files will be created, design_file.txt and design_file.txt.sum. The design_file.txt is the design file for job_expression.pl, the design_file.txt.sum is the design file for summarize4expression.py. To get more detialed help, please run: python3 design4expression.py -h 
'''


import csv
import os


# Import tools.py and rcctools.py in different cases
import pathlib

rootDir = os.path.dirname(__file__)

fullPath2tools = os.path.join(rootDir, "tools.py")
if os.path.isfile(fullPath2tools):
	import tools
elif os.path.isfile(os.path.join(rootDir, "script", "tools.py")):
	# assume tools.py in script subdirectory
	from script import tools

fullPath2rcctools = os.path.join(rootDir, "rcctools.py")
if os.path.isfile(fullPath2rcctools):
	import rcctools
elif os.path.isfile(os.path.join(rootDir, "script", "rcctools.py")):
	# assume rcctools.py in script subdirectory
	from script import rcctools


# Refer to "/project/shared/xiao_wang/software/rnaseqDE/job_expression.pl" for the details about format.
columnNames4sum = ('Patient_ID', 'Tumor_sample_ID', 'Path2output')

def mainFunc(args4mainFunc):
	manifest_file = args4mainFunc['manifest_file']
	rows4sample = tools.readDictCsvFile(manifest_file, delimiter='\t')
	
	patientDic = {}
	for row in rows4sample:
		Patient_ID = row['Patient_ID'].strip()
		Root_tree = row['Root_tree'].strip()
		Data_type = row['Data_type'].strip()
		if Patient_ID not in patientDic.keys():
			patientDic[Patient_ID] = {}
		sampleID = (Root_tree, Data_type)
		patientDic[Patient_ID][sampleID] = row

	result_path = args4mainFunc['result_path']
	designList = []
	designList4sum = [columnNames4sum]
	for Patient_ID,sampleDic in patientDic.items():
		for sampleID, sample in sampleDic.items():
			Root_tree, Data_type = sampleID
			if Data_type not in ('RNA', 'RNAseq', 'RNA-seq'):
				continue
			Path2output = os.path.join(result_path, Root_tree)
			Seq_ID = sample['Seq_ID'].strip()
			Path = sample['Path'].strip()
			if 'Species' in sample.keys() and sample['Species'] is not None:
				Species = sample['Species'].strip()
			else:
				Species = ''

			if not os.path.isdir(Path):
				print('Warning, directory ({}) not found'.format(Path))
				print('Skip sample', Root_tree)
				continue
			# search for files in directory Path, where the basename of file contains string Seq_ID and ends with end
			if '.fastq' in Seq_ID or '.gz' in Seq_ID:
				end = ''
			else:
				end = tools.fastqends()
			fastqfiles = tools.getfiles(Path, Seq_ID, endswith=end)
			fastqfiles.sort()
	
			nfiles = len(fastqfiles)
			if nfiles == 2:
				f1, f2 = fastqfiles
				seqfile1 = os.path.join(Path, f1)
				seqfile2 = os.path.join(Path, f2)
			elif nfiles == 1:
				f1 = fastqfiles[0]
				seqfile1 = os.path.join(Path, f1)
				seqfile2 = ''
			else:
				print('Warning: {} (!=2) files containing Seq_ID {} in {}: {}'.format(
					nfiles, Seq_ID, Path, fastqfiles))
				print('Ignore the sample', Root_tree)
				continue

			pdx = rcctools.pdx(sample)
			row = (seqfile1, seqfile2, Path2output, pdx)
			designList.append(row)
			share_path = args4mainFunc['share_path']
			if (share_path == ''):
				row4sum = (Patient_ID, Root_tree, Path2output)
			else:
				Path2output = os.path.join(share_path, Root_tree)
				row4sum = (Patient_ID, Root_tree, Path2output)
			designList4sum.append(row4sum)
	#designList = tools.removeExistingDir(designList, lambda x: x[-2])

	designFile = args4mainFunc['design_file']
	if len(designList) > 1:
		tools.makeDirTree(result_path)
		tools.writeCsvFile(designFile, designList, delimiter='\t')
	else:
		print('No new RNA sample is found, no new design file ({}) is created'.format(designFile))
	designFile4sum = designFile + '.sum'
	if len(designList4sum) > 1:
		tools.writeCsvFile(designFile4sum, designList4sum, delimiter='\t')
	else:
		print('No RNA sample is found, no design file for summarization ({}) is created'.format(designFile4sum))
	

if __name__ == "__main__":
	import argparse
	import textwrap

	# Parse command line arguments
	descriptStr = '''The script takes a manifest file as the input, and creates two design files for job_expression.pl and summarize4expression.py.'''
	parser = argparse.ArgumentParser(prog='design4expression.py', description = textwrap.dedent(descriptStr), 
			formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument(
		'manifest_file', 
		type=str,
		default='',
		help = 'Input, manifest file of samples, must have header line with column names, see example_data/manifest.txt',
		)
	parser.add_argument(
		'result_path', 
		type=str,
		default='',
		help = 'Root directory in which the sample folders containing the output of expression.pl are placed',
		)
	parser.add_argument(
		'design_file', 
		type=str,
		default='',
		help = 'Output, the design file for job_expression.py. Another design file with suffix .sum for summarize4expression.py will be created as well',
		)
	parser.add_argument(
		'--share_path', 
		required = False,
		type=str,
		default='',
		help = 'Root path to share folder where the results are shared with others. If provided, it replaces result_path in the design file for summarize4expression.py',
		)
	

	args = parser.parse_args()

	args4mainFunc =	{
				'manifest_file': args.manifest_file.strip(),
				'result_path': args.result_path.strip(),
				'design_file': args.design_file.strip(),
				'share_path': args.share_path.strip(),
			}

	mainFunc(args4mainFunc)

