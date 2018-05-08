import sys, argparse, tom, multiprocessing, os, pickle, random, glob, gzip
import tom_slim as ts
import site_frequency_spectrum as sfs_tools
from collections import Counter
from tom import brace

def orgPolyDict(organ_mutations,N):
	sfsDict = {}
	for i in organ_mutations.keys():
		orgDict = {}
		orgDictRaw = {}
		for q in organ_mutations[i]:
			if q[1] not in orgDictRaw.keys():
				orgDictRaw[q[1]] = [int(q[7])]
			else:
				orgDictRaw[q[1]].append(int(q[7]))
		for o in orgDictRaw.keys():
			orgDict[o] =sfs_tools.SFS_from_all_frequencies(orgDictRaw[o],N) 
#			print orgDict[o]
		sfsDict[i] = orgDict
	return sfsDict

def orgFixDict(fixed_organs):
	fixD = {}
	for i in fixed_organs.keys():
		orgFixD = Counter()
		for o in fixed_organs[i]:
			orgFixD[o[1]]+=1
		fixD[i] = orgFixD
	return fixD
			
def combinePolyFix(poly,fix):
#	print poly
#	print fix.keys(),poly.keys()
	for i in fix.keys():
#		print i
		if i not in poly.keys(): continue
		if len(poly[i]) == 0: continue
		for  j in fix[i]:
#			print j
			try:
				if j not in poly[i].keys():
					
					poly[i][j] = [0]*len(poly[i][poly[i].keys()[0]])
					poly[i][j][-1] += fix[i][j]
				else:
					poly[i][j][-1] += fix[i][j]
			except KeyError:continue		
	return poly

def combineElements(polyfix,lengths):
	elDict = {}
	for i in polyfix.keys():
		sfs = []
		try:
			for j in polyfix[i].keys():
				if i == 'g1' and j == 'm12': ## for the case of exons
					elDict['syn'] = polyfix[i][j]
					continue
				if len(sfs) == 0:
					sfs = polyfix[i][j]
				else:
					sfs = sfs_tools.merge_SFS(sfs,polyfix[i][j])
		except KeyError:continue
		elDict[i] = sfs

	for i in elDict.keys():
#		print i
		zero = lengths[i] - sum(elDict[i])
		try:
			elDict[i][0] = zero
		except IndexError:
			return elDict
	return elDict


def parseLengths(lengthDictRaw):
	lengthDict = {}
	for i in lengthDictRaw.keys():
		if i == 'g1':
			raw = lengthDictRaw[i]
			syn = raw * 0.25
			nonsyn = raw * 0.75
			lengthDict['syn'] = syn
			lengthDict['g1'] = nonsyn
		else:
			lengthDict[i] = lengthDictRaw[i]
	return lengthDict

def get_sfs_dict(slim_input,num=-1):
	data = [i.strip() for i in gzip.open(slim_input).readlines()]
	x = ts.slim(data,fixed = True,give_genomes=True)
#	print x.genomes
	if not x.sanity:
		return [None,None]
	thresh = x.N*10
	
	polyDict = orgPolyDict(x.organ_mutations(),x.sampleN)
	lengthDict = parseLengths(x.organ_lengths())
#	print lengthDict
#	print polyDict
	fixedDict = x.organ_fixed(threshold = int(x.N)*20)
	
	fixD = orgFixDict(fixedDict)
	
	polyfix = combinePolyFix(polyDict,fixD)
#	print lengthDict
	elDict = combineElements(polyfix,lengthDict)	
	print 'processed ' + x.name
	return [x.name,elDict]

def get_sfs_dict_from_sample(slim_input):
	data = [i.strip() for i in gzip.open(slim_input).readlines()]
	x = ts.slim(data,fixed = True,give_genomes=True)
	if not x.sanity:
		return [None,None]
#	print x.name
	genomes = x.genome_dict()
	mutations = x.mutations_dict()
	lengthDict = parseLengths(x.organ_lengths())
	individuals = [random.choice(genomes.keys()) for i in range(20)]
#	print individuals
#	if x.name == 	'/exports/csce/eddie/biology/groups/eddie_biology_ieb_keightley/toms_simulations/updated_DFE/longRuns/full_usfs/configs/3381.temp.slim':
#		individuals = ['p1:1398', 'p1:1646', 'p1:297', 'p1:165', 'p1:999', 'p1:1451', 'p1:982', 'p1:973', 'p1:615', 'p1:832', 'p1:12', 'p1:1109', 'p1:1137', 'p1:496', 'p1:164', 'p1:412', 'p1:1687', 'p1:1373', 'p1:72', 'p1:39']
	muts_by_organ = x.organ_mutations()
	new_muts = Counter()
	for g in individuals:
		for m in genomes[g]:
			new_muts[m] +=1
	
	polyDict = {} 
	for h in muts_by_organ.keys():
		mTypeDict = {}
		for m in muts_by_organ[h]:
			if new_muts[m[0]] == 0: continue
			if m[1] not in mTypeDict.keys():
				mTypeDict[m[1]] = [new_muts[m[0]]]
			else:
				mTypeDict[m[1]].append(new_muts[m[0]])
#		print h, mTypeDict
		mPoly = {}
		for k in mTypeDict.keys():
			mPoly[k] = sfs_tools.SFS_from_all_frequencies(mTypeDict[k],20) 
		polyDict[h] = mPoly

#	print '!', x.name
	fixedDict = x.organ_fixed(threshold = int(x.N)*20)

	fixD = orgFixDict(fixedDict)
	polyfix = combinePolyFix(polyDict,fixD)
	elDict = combineElements(polyfix,lengthDict)	

	
	print 'processed ' + x.name
	return [x.name,elDict]



def main():
	parser = argparse.ArgumentParser(description="Get the SFS for each mutation type in the simulation")
	parser.add_argument("-i","--input", 
		required = True,
		dest = "input", 
		type =str, 
		help = "Give the name of the gzipped, SLiM output")

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output", 
		type =str, 
		help = "the name of the output, this script will gzip that output automatically")
	parser.add_argument("-t","--procs", 
		required = True,
		dest = "procs", 
		type =int, 
		help = "The number of processors to use")
	parser.add_argument("-n","--num", 
		required = False,
		dest = "num", 
		type =int, 
		help = "The number of individuals to downsample to",
		default = -1)

#	slims =['1639.out.gz','516.out.gz','3791.out.gz','2519.out.gz','3277.out.gz','2615.out.gz','3649.out.gz','3742.out.gz','11.out.gz','128.out.gz','1967.out.gz','325.out.gz','1480.out.gz','1271.out.gz','267.out.gz','3935.out.gz','2922.out.gz','1134.out.gz']
	slims =['516.out.gz','3791.out.gz','2519.out.gz','3277.out.gz','128.out.gz','1967.out.gz','325.out.gz','1271.out.gz','267.out.gz','1134.out.gz']

	slims2 = ['Ne-anc/runs/'+ll for  ll in slims if ll == '1134.out.gz']
	args = parser.parse_args()
	finalDict = {}
	count = 0
	if args.num > 0:
		sfsDictFunc = get_sfs_dict_from_sample
	elif args.num == -1:
		sfsDictFunc = get_sfs_dict
	

	if args.procs == 1:	
		print 'Reading all SLiM output'
		for i in glob.glob(args.input+'/*'):
#		for i in ts.slim_reader_gzip(args.input):
			count +=1
			print i
		
			runDict = sfsDictFunc(i)
			
			finalDict[runDict[0]]=runDict[1]
			
	elif args.procs > 1:
		p = multiprocessing.Pool(args.procs)
		tempSlims = [g for g in glob.glob(args.input+'/*') if g in slims2]

		print 'Reading all SLiM output'
#		results = p.map( sfsDictFunc, tempSlims)
		results = p.map( sfsDictFunc, glob.glob(args.input+'/*'))
##		results = p.map( sfsDictFunc, ts.slim_reader_gzip(args.input))
		for j in results:
			if j[0] == None:continue 
			finalDict[j[0]]=j[1]
#	print '!'	
	pickle_jar = open(args.output,"wb")
#	print '!!'	
	pickle.dump( finalDict , pickle_jar )
        pickle_jar.close()
	os.system("gzip "+args.output)

if '__name__':
	main()
