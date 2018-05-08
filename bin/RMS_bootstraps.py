### Use bootstraps to get confidence intervals around the RMS of the trough in diversity 

import pickle, argparse, gzip
import site_frequency_spectrum as SFS

def main():
	parser = argparse.ArgumentParser(description="Takes the pickle file of intervals around functional elements and bootstraps them, gets the pattern of diversity, then calculates RMS from an analysis file")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the pickle that you want to analyse")
	parser.add_argument("-c","--comparison", 
		required = True,
		dest = "comparison",
		type = str, 
		help = "The analysis files that you want to compare the simulations to")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")

	args = parser.parse_args()

	data = pickle.load(gzip.open(args.input))

	print data.keys()


if '__name__':
	main()


