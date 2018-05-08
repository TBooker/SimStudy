import argparse, pandas as pd, math
from tom import brace

def compare_slice(sliced, comparison, element):

	comparison['pi_div_t2']  = comparison['pi'] / comparison['t2'] 
	if element == 'exon':
		dist_for_comp = 75000
		comp1 = comparison['pi'][comparison["dist"].abs() > dist_for_comp].mean() 
		comparison['pi_red']  = comparison['pi'] / comp1

		comp2 = sliced['pi'][sliced["mid"].abs() > dist_for_comp].mean() 
		sliced['pi_red']  = sliced['pi'] / comp2

		comparison = comparison.set_index(comparison['dist'].values).sort_index()
		sliced = sliced.set_index(sliced['start'].values -1).sort_index()


	elif element == 'cne':
		dist_for_comp = 4000
		mean_t2 = comparison['pi'][comparison["dist"].abs() > dist_for_comp].mean() 
		comp1 = comparison['pi_div_t2'][comparison["dist"].abs() > dist_for_comp].mean() 
		comparison['pi_red']  = comparison['pi_div_t2'] / comp1

		comp2 = sliced['pi'][sliced["mid"].abs() > dist_for_comp].mean() 
		sliced['pi_red']  = sliced['pi'] / comp2

		comparison = comparison.set_index(comparison['dist'].values).sort_index()
		sliced = sliced.set_index(sliced['mid'].values).sort_index()

#	print sliced
#	print comparison	
	return math.sqrt(((comparison['pi_red'] - sliced['pi_red']) **2).mean())

def main():

	parser = argparse.ArgumentParser(description="From the input file, calculates the ")

	parser.add_argument("-i","--input", 
			required = True,
			dest = "input",
			type =str, 
			help = "The name of the file that contains the SLiM output")
	parser.add_argument("-o","--output", 
			required = True, 
			metavar ="output", 
			type = str, 
			help = "What name do you want to gvve to the output file full of stats?")
	parser.add_argument("-c","--comparison", 
			metavar ="comparison", 
			required = True,
			type = str, 
			help = "Give the name of the file with the observed data in it, will use to get RMS")
	parser.add_argument("--cne", 
			required = False,
			action = 'store_true', 
			help = "Are you comparing CNE? If so, need to do a little transform")

	args = parser.parse_args()

	boots = pd.read_csv(args.input)	
	comparison = pd.read_csv(args.comparison)
	number_of_replicates = len(set(boots['label']))
	brace()
	if args.cne:
		element = 'cne'
	else:
		element = 'exon'
	print 'analysing:',element
	rms = []

	for i in set( boots['label'] ):
		sliced = boots[boots['label'] == i].copy()
		rms_for_slice = compare_slice(sliced, comparison, element)
		rms.append( rms_for_slice )
		

	print 'rms_lower:', sorted(rms)[int(number_of_replicates*0.025)]
	print 'rms_median:', sorted(rms)[int(number_of_replicates*0.5)]
	print 'rms_upper:', sorted(rms)[int(number_of_replicates*0.975)]


if '__name__':
	main()
