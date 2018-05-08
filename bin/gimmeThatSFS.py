import sys, argparse, tom, gzip, random
import cPickle as pickle
import site_frequency_spectrum as SFS_tools






def main():
	parser = argparse.ArgumentParser(description="Take the dict of SFS for different SLiM runs and make a composite SFS or bootstrapped SFS")
	parser.add_argument("-i","--input", 
		required = True,
		dest = "input", 
		type =str, 
		help = "Give the name of the gzipped pickle file")

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output", 
		type =str, 
		help = "the name of the output")
	parser.add_argument("-b","--boots", 
		required = False,
		action = 'store_true',
		help = "Do you want a bootstrapped SFS?",
		default = False)

	args = parser.parse_args()
	
	pickle_jar = pickle.load(gzip.open(args.input,"rb"))
	

	names = [p for p in pickle_jar.keys() if p]
	for i in names:
		print i
	if args.boots:
		boots = [names[random.randint(0,len(names)-1)] for i in range(len(names))]
		names = boots

	bigDict = {}
	for n in names:
		z = pickle_jar[n]
		for j in z.keys():
                        if not z[j]: continue
			if j not in bigDict.keys():
                                bigDict[j] = z[j]
			else:
				bigDict[j] =SFS_tools.merge_SFS(z[j],bigDict[j]) 
	output = open(args.output,'w')
	for i in bigDict.keys():
		sfs = map(str,bigDict[i])
		output.write(i+'\n'+' '.join(sfs)+'\n')
	output.close()

if "__name__":
	main()
