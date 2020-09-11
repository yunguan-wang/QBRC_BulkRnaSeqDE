import os
import io
import csv
import statistics
import shutil

if __package__ is None or __package__ == '':
	import constants
else:
	from . import constants


# Return a copy of the iterable where strip() is applied on each string elelment with the leading 
# and trailing characters if the element is not a iterable, else cleanrow() will be applied to the 
# element recursively.
# 
# row: iterable, e.g. ['abc', 'printer ', ' hello world ']
# 
# Return newrow.
# newrow: cleaned row, every element of which has leading and trailing characters removed, the 
#	characters are specified by 'chars'.
def cleanrow(row, chars):
	newrow = []
	for v in row:
		if isinstance(v, str):
			newv = v.strip(chars)
		elif len(v) > 1:
			newv = cleanrow(v, chars)
		newrow.append(newv)
	return newrow


# read in all rows in csv file, including the empty rows
def readCsvFileAll(csvfile, delimiter='\t'):
	rows = []
	with open(csvfile, newline='') as fp:
		reader = csv.reader(fp, delimiter=delimiter)
		for row in reader:
			if len(row) == 0:
				row = ['']
			rows.append(row)
	return rows

# read in all rows except the empty rows in csv file
def readCsvFile(csvfile, delimiter='\t'):
	rows = []
	with open(csvfile, newline='') as fp:
		reader = csv.reader(fp, delimiter=delimiter)
		for row in reader:
			# skip the rows without values
			rowlength = [len(x.strip()) for x in row if x != None]
			if len(rowlength) < 1 or max(rowlength) < 1:
				continue
			rows.append(row)
	return rows

# Skip the first line 
def readCsvFileSkipLine(csvfile, delimiter='\t'):
	rows = []
	with open(csvfile, newline='') as fp:
		reader = csv.reader(fp, delimiter=delimiter)
		next(reader)
		for row in reader:
			# skip the rows without values
			rowlength = [len(x.strip()) for x in row if x != None]
			if len(rowlength) < 1 or max(rowlength) < 1:
				continue
			rows.append(row)
	return rows

# Read each row in .csv with head line into a dictionary with the column name as the key of the correponding element
# in the dictionary.
#
# Input
# csvfile:	character string, full-path file name pointing to a csv file
# delimiter:	character, delimiter of columns in csv file
#
# Return dicList
# dicList:	[dic, ..., dic]
#  dic:		ordereddic{key:value, ..., key:value}
#
def readDictCsvFile(csvfile, delimiter=',', fieldnames=None):
	dicList = []
	with open(csvfile, newline='') as fp:
		reader = csv.DictReader(fp, delimiter=delimiter, fieldnames=fieldnames)
		for row in reader:
			# skip the rows without values
			if max([len(x.strip()) for k,x in row.items() if k != None and x != None]) < 1:
				continue
			dicList.append(row)
	return dicList


# Input
# csvfile:	character string, full-path file name pointing to a csv file
# fieldnames:	[filedname, ...], e.g. ['ID', 'value1', 'value2', ...]
# dicList:	[dic, ..., dic]
#  dic:		ordereddic{key:value, ..., key:value}
# 
def writeDictCsvFile(csvfile, fieldnames, dicList, delimiter=','):
	with open(csvfile, 'w', newline='') as fp:
		#writer = csv.DictWriter(fp, fieldnames=fieldnames, extrasaction='ignore', delimiter=delimiter)
		writer = csv.DictWriter(fp, fieldnames=fieldnames, lineterminator='\n', delimiter=delimiter)
		writer.writeheader()
		writer.writerows(dicList)

def writeDictCsvFileWithoutHeader(csvfile, fieldnames, dicList, delimiter=','):
	with open(csvfile, 'w', newline='') as fp:
		writer = csv.DictWriter(fp, fieldnames=fieldnames, lineterminator='\n', delimiter=delimiter)
		#writer.writeheader()
		writer.writerows(dicList)

# Input
# csvfile:	character string, full-path file name pointing to a csv file
# rowList:	[row, ...]
#  row:		(item, ...) or [item, ...]
def writeCsvFile(csvfile, rowList, delimiter=','):
	with open(csvfile, 'w', newline='') as fp:
		spamwriter = csv.writer(fp, delimiter=delimiter, lineterminator='\n')
		spamwriter.writerows(rowList)

def intersection(p1, p2):
	a, b = p1
	c, d = p2
	# Must be satisfied: a <= b and c <=d 
	overlap = min(b, d) - max(a, c) + 1
	return overlap

# geneDic[gene][sample]='' for sample not in geneDic
def addSamples(geneDic, samples):
	for gene, gene2samples in geneDic.items():
		for sample in samples:
			if sample not in gene2samples.keys():
				geneDic[gene][sample] = ''
	
# output gene-sample matrix
# geneDic: {gene:gene2samples, ...}
# gene2samples: {sample:value, ...}
def outputSumFile(geneDic, sumFile, samples):
	row4title = [''] + sorted(samples)
	#str4row4title = '\t'.join(row4title)
	rows = [row4title]
	for gene, gene2samples in geneDic.items():
		row = [gene]
		for sample in sorted(gene2samples.keys()):
			counts = gene2samples[sample]
			row.append(counts)
		rows.append(row)
	writeCsvFile(sumFile, rows, delimiter='\t')

# get the list of chromosome arms from armLoc
def getArms(armLoc):
	armids = []
	for chr,arms in armLoc.items():
		for arm in arms.keys():
			armids.append(arm)
	return armids

# Return armLoc
# armLoc: {chr:arms, ...}
# chr: chromosome ID
# arms:{arm:(start, end, bands), ...}
# arm: arm ID, like chr1_p, chr1_q
# bands: [band, ...]
# band: (start, end, bandID)
# bandID: band ID, e.g. p36.33 and q11
def getArmLoc(cytoBandFile):
	armLoc = {}
	with open(cytoBandFile, 'r') as fp:
		rawArms= readDictCsvFile(cytoBandFile, delimiter='\t')
	for row in rawArms:
		# remove the band without arm information
		if len(row['band'].strip()) == 0:
			continue
		chr = row['chromosome'].strip()
		bandID = row['band'].strip()
		start = int(row['start'].strip())
		end = int(row['end'].strip())
		if chr not in armLoc:
			armLoc[chr] = {}
			armName = bandID[0]
			start4arm = start
			bands = []
		if bandID[0] != armName:
			armName = bandID[0]
			start4arm = start
			bands = []
		end4arm = end
		band = (start, end, bandID)
		bands.append(band)
		armid = chr + '_' + armName
		armLoc[chr][armid] = (start4arm, end4arm, bands)
	return armLoc

# Return armLoc
# armLoc: {chr:arms, ...}
# chr: chromosome ID
# arms:{arm:(start, end, bands), ...}
# arm: arm ID with band info, like chr1_p36.33, chr1_q11
# bands: [band, ...]
# band: (start, end, bandID)
# bandID: band ID, e.g. p36.33 and q11
def getArmLocByBand(cytoBandFile):
	armLoc = {}
	with open(cytoBandFile, 'r') as fp:
		rawArms= readDictCsvFile(cytoBandFile, delimiter='\t')
	for row in rawArms:
		# remove the band without arm information
		if len(row['band'].strip()) == 0:
			continue
		chr = row['chromosome'].strip()
		bandID = row['band'].strip()
		start = int(row['start'].strip())
		end = int(row['end'].strip())
		if chr not in armLoc:
			armLoc[chr] = {}
			armName = bandID
			start4arm = start
			bands = []
		if bandID != armName:
			armName = bandID
			start4arm = start
			bands = []
		end4arm = end
		band = (start, end, bandID)
		bands.append(band)
		armid = chr + '_' + armName
		armLoc[chr][armid] = (start4arm, end4arm, bands)
	return armLoc

# get the list of IDs in dictionary
# armLoc: {chr:arms, ...}
# arms: {arm:(start, end, bands), ...}
# arm: arm ID (e.g. chr1_p, chr1_q) or band ID (e.g. chr1_p36.33 and chr1_q11)
def getArms(armLoc):
	armids = []
	for chr,arms in armLoc.items():
		for arm in arms.keys():
			armids.append(arm)
	return armids


# geneLoc: {gene:loc, ...}
# loc: (chr, start, end)
def getGeneLoc(rows4refFlat):
	geneLoc = {}
	for row in rows4refFlat:
		gene = row['name2']
		chr = row['chr']
		start = int(row['txStart'].strip())
		end = int(row['txEnd'].strip())
		geneLoc[gene] = (chr, start, end)
	return geneLoc

# Return armDic
# armDic: {arm:arm2samples, ...]
# arm2samples: {sample:ncnv, ...}
# ncnv: float, average of the cnv where ncnv is the normalized cnv of genes in a sample
def summarizeByArm(geneDic, geneLoc, armLoc, arms, samples):
	armDic = {}
	# initailize armDic
	for arm in arms:
		armDic[arm] = {}
		for sample in samples:
			armDic[arm][sample] = []
	# fill events in armDic
	arm2genes2samples = {}
	for gene,gene2samples in geneDic.items():
		if gene not in geneLoc.keys():
			continue
		chr, start, end = geneLoc[gene]
		if chr not in armLoc.keys():
			continue

		segmentbd = (start, end)
		segmentLen = max(segmentbd) - min(segmentbd) + 1

		for arm,arminfo in armLoc[chr].items():
			armbd = arminfo[:2] 
			overlap = intersection(segmentbd, armbd)
			armLen = max(armbd) - min(armbd) + 1
			# we say the segment is located in the arm region if overlapBetweenSegmentAndArm/min(segmentLen,armLen) > 0.5
			if overlap / min(segmentLen, armLen) > 0.5:
				if arm not in arm2genes2samples.keys():
					arm2genes2samples[arm] = []
				arm2genes2samples[arm].append([gene2samples, segmentbd])
	for arm,genes2samples in arm2genes2samples.items():
		for gene2samples2segBd in genes2samples:
			gene2samples, segmentbd = gene2samples2segBd
			for sample,cnv in gene2samples.items():
				seg = (cnv, segmentbd)
				armDic[arm][sample].append(seg)
			
	# Averaging cnv values of all CNV records of the specific arm occuring in the specific sample
	# because there might be more than one CNV segment located in the specific arm.
	# Note: the value of cnv of the specific arm occuring in the specific tumor sample is set to '' if n == 0
	for arm,arm2samples in armDic.items():
		for sample,segs in arm2samples.items():
			# no CNV events occurs in the specific chromosome arm in the specific sample
			n = len(arm2samples[sample])
			if n == 0:
				armDic[arm][sample] = ''
			# at least one CNV events occurs in the specific chromosome arm in the specific sample
			elif n >= 1:
				# weighted average by segment length
				sum4cnvWeightedByLen = sum([
							seg[0] * (seg[1][1]-seg[1][0]+1) 
							for seg in segs
							])
				nWeightedByLen = sum([
							seg[1][1]-seg[1][0]+1 
							for seg in segs
							])
				armDic[arm][sample] = sum4cnvWeightedByLen / nWeightedByLen

	"""
	# normalize cnv in each sample
	cnv4sample = {}
	for arm,arm2samples in armDic.items():
		for sample,cnv in arm2samples.items():
			if sample not in cnv4sample:
				cnv4sample[sample] = []
			if str(cnv) != '':
				cnv4sample[sample].append(cnv)
	# m: {sample:median, ...}
	# median: median of cnv of all CNV segments occuring in a sample
	m = {}
	for sample,cnvs in cnv4sample.items():
		m[sample] = statistics.median(cnvs)
	# normalize the cnv
	for arm,arm2samples in armDic.items():
		for sample,cnv in arm2samples.items():
			if str(cnv) == '':
				armDic[arm][sample] = m[sample]
			else:
				armDic[arm][sample] = cnv - m[sample]
	"""
	return armDic

# remove the row if all values in a row are empty string or blank space.
def removeBlankRow(armDic):
	armDicNew = {}
	for arm,arm2samples in armDic.items():
		values = set(str(v).strip() for v in arm2samples.values())
		if len(values) == 1 and values.pop() == '':
			continue
		armDicNew[arm] = arm2samples
	return armDicNew

# list possible suffice of file names for fastq files
def fastqends():
	#return ['.gz']
	#return ['.fastq.gz']
	ends = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']
	return ends

def bamends():
	return ['.bam']

# list possible pattern of file names which contain Seq_ID
def getSeqids(Seq_ID):
	return [Seq_ID]

	'''
	ids = []
	seps = ['.', '_', '-']
	rs = ['r', 'R']
	pairs = ['1', '2']
	for sep in seps:
		for r in rs:
			for number in pairs:
				id = Seq_ID + sep + r + number
				ids.append(id)
	return ids
	'''

# list possible patterns of character strings which contain seeds
def getPatterns(seeds = ''):
	ids = []
	if type(seeds) is str:
		seeds = [seeds]
	seps = ['.', '_', '-']
	for seed in seeds:
		for sep in seps:
			id = sep + seed + sep
			ids.append(id)
	return ids

# search for files in directory path, where the basename of the file contains str and ends with endswith
# Return files
# files: [f, ...], f is the base file name where the path to file is removed
# id: character string or a list of character strings
# endswith: character string or a list of character strings
# Note: str.endswith('') always return True, namely, getfiles() does not check the end of the file name
def getfiles(path, ids, endswith=''):
	seps = ('.', '_', '-', 'SAM')
	files = []
	if type(ids) is str:
		ids = [ids]
	if type(endswith) is str:
		endswith = [endswith]
	ids = tuple(ids)
	endswith = tuple(endswith)
	with os.scandir(path) as items:
		for item in items:
			if item.is_file():
				name = item.name
				if not name.startswith('.'):
					for id in ids:
						if id in name:
							items = name.split(id)
							if ((len(items[0]) == 0 or items[0].endswith(seps)) and
								(len(items[1]) == 0 or items[1].startswith(seps))):
								idflag = True
								break
							# dirty code to retrieve file names with id
							if ((len(items[0]) == 0 or items[0].endswith(seps)) and
								(len(items[1]) == 0 or items[1][0] == 'R')):
								idflag = True
								break
							# dirty code to retrieve file names with id
					else:
						idflag = False
					#for end in endswith:
					#	if name.endswith(end):
					#		endflag = True
					#		break
					if name.endswith(endswith):
						endflag = True
					else:
						endflag = False
					if idflag == True and endflag == True:
						# temporary and dirty codes to remove unwanted files
						#if not name.startswith(ids):
						#	continue
						if '_I1_' in name:
							continue
						if '_val_' in name:
							continue
						# temporary and dirty codes to remove unwanted files
						files.append(name)
	return files

# check if any valid NGS sequence file (fastq) with size > 0 exist at the specified path
def fastq_seq_files_exist(Seq_ID, Path):
	Seq_ID = getSeqids(Seq_ID)
	end = fastqends()
	files = getfiles(Path, Seq_ID, endswith=end)
	files.sort()
	nfiles = len(files)
	if nfiles < 1:
		flag = False
	else:
		if nfiles > 2:
			print('Warning: {} (>2) fastq files with {} at {}'.format(
				nfiles, Seq_ID, Path))
		for f in files:
			f = os.path.join(Path, f)
			if os.path.isfile(f) and os.stat(f).st_size > 0:
				fsize = os.stat(f).st_size
				if fsize < constants.GB:
					print("Warning: small ({}MB) fastq file {} at {}".format(
						fsize/constants.MB, f, Path))
				flag = True
				break
		else:
			flag = False
	if flag:
		return True
	else:
		return False

# check if any valid NGS sequence file (either fastq or bam file) with size > 0 exist at the specified path
def seq_files_exist(Seq_ID, Path):
	Seq_ID = getSeqids(Seq_ID)
	end = fastqends() + bamends()
	files = getfiles(Path, Seq_ID, endswith=end)
	files.sort()
	nfiles = len(files)
	if nfiles < 1:
		flag = False
	else:
		for f in files:
			f = os.path.join(Path, f)
			if os.path.isfile(f) and os.stat(f).st_size > 0:
				flag = True
				break
		else:
			flag = False
	if flag:
		return True
	else:
		return False


# remove from list 'files' the items containing substrings in patterns, return filesRefined
# filesRefined: [str, ...]
#
# files: [str, ...]
def removeitems(files, patterns=[]):
	filesRefined = []
	for f in files:
		for p in patterns:
			if p in f:
				break
		else:
			filesRefined.append(f)
	return filesRefined

# same as getfiles() but check and return directories instead of files
def getdirs(path, str, endswith=''):
	dirs = []
	with os.scandir(path) as items:
		for item in items:
			if item.is_dir():
				name = item.name
				if not name.startswith('.') and name.endswith(endswith) and str in name:
					dirs.append(name)
	return dirs

# Return True if it is a normal sample, else False
# Assume that:
#	1) Root_tree of normal sample should be patient_ID + 'N' + optional non-alpha character, e.g. AS10N, XP66N2
#	#2) Source contains word 'normal', e.g. 'kidney normal'
def isnormalsample_v1(sample):
	Root_tree = sample['Root_tree'].strip()
	#Source = sample['Source'].strip()
	Patient_ID = sample['Patient_ID'].strip()
	# Source contains word 'normal', e.g. 'kidney normal'
	#if 'normal' in Source.lower():
	# Root_tree of normal sample should be patient_ID + 'N' + optional non-alpha character, e.g. AS10N, XP66N2
	Root_treeStart4normal = Patient_ID + 'N'
	len4Root_treeStart4normal = len(Root_treeStart4normal)
	if Root_treeStart4normal == Root_tree[:len4Root_treeStart4normal]:
		return True
	else:
		return False

# Return True if it is a normal sample, else False
# Assume that:
#	1) The tail of Root_tree of normal sample is 'N' + optional character
#	   e.g. AS10N for patient AS10, XP66N2 for patient XP66, 005N for patient 92, XP409Na for patient XP409
#	AND
#	2) sample['Source'] contains word 'normal', e.g. 'kidney normal'
def isnormalsample(sample):
	Root_tree = sample['Root_tree'].strip()
	Source = sample['Source'].strip().lower()
	if 'N' in Root_tree and 'normal' in Source:
		return True
	else:
		Patient_ID = sample['Patient_ID'].strip()
		if 'N' in Root_tree or 'n' in Root_tree:
			print("Warning: 'N' or 'n' in Root_tree of tumor sample {} for patient {}".format(Root_tree, Patient_ID))
		elif 'normal' in Source:
			print("Warning: 'Source' {} of tumor sample {} for patient {} contains normal".format(Source, Root_tree, Patient_ID))
		return False
	

# Return listRefined
# listRefined: refined designList with some rows removed
# designList: [row, ...]
# row: diectionary or iterable such as list or tuple
def removeExistingDir(designList, func4path):
	listRefined = []
	for row in designList:
		path = func4path(row)
		if os.path.isdir(path):
			print('Ignore existing', path)
			continue
		listRefined.append(row)
	return listRefined


# Check if the file exist or not
# Return True if existing, else False
# f: full path to a file
def checkFileExist(f):
	if os.path.isfile(f):
		return True
	else:
		return False

# transpose a matrix
def transpose_matrix(matrix):
    matrix = list(map(list, zip(*matrix))) # assign result 
    return matrix # return transposed matrix

# makes all intermediate-level directories needed to contain the leaf directory
def makeDirTree(path):
	if not os.path.exists(path):
		try:
			os.makedirs(path)
		except OSError as error:
			if error.errno != errno.EEXIST:
				raise
	elif os.path.isfile(path):
		e = "Fail to create the directory {}, it exists and is a regular file.".format(path)
		raise RuntimeError(e)

# remove non-empty directory
def removeDirTree(path):
	try:
		shutil.rmtree(path)
	except OSError as e:
		print("Error: %s : %s" % (path, e.strerror))

