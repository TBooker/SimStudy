### Script to get use the Kimura equation for the fixation of a selected allele
import argparse, math, sys
import pandas as pd
#double kimura_fixation_prob1(double s, double N)
#{
#   double num, denom;
#   if (s == 0.0)
#   {
#      return 1.0/(2.0*N);
#   }
#   else if ((s < 0.0)&&(N*s > -0.0000001))  // Trap numerical problems
#   {
#      num = 1.0;
#      denom = 2.0*N;
#   }
#   else if ((s > 0.0)&&(N*s < 0.0000001))
#   {
#      num = 1.0;
#      denom = 2.0*N;
#   }
#   else
#   {
#      num = 1 - exp(-s);
#      denom = 1 - exp(-2*N*s);
#   }
#//   printf("s %14.14lf N %lf num %lf denom %lf num/denom %lf\n",  s, N, num, denom, num/denom);
#//   monitorinput();
#   return num/denom;
#}

def kimura_fixation_prob(s,N):
	if s == 0.0:
		return 1.0/(2*N)
	elif s < 0.0 and N*s > -0.00000000001:
		num = 1.0
		denom = 2.0*N
	elif s > 0.0 and N*s < 0.00000000001:
		num = 1.0
		denom = 2.0*N
	else:
		num = 1 - math.exp(-s)
		denom = 1 - math.exp(-2*N*s)
	return num/denom
def evaluateKaKs(pa,sa,dfe, classes):
	Nw = dfe['Nw']
	try:
		Na = dfe['N_a']
	except KeyError:
		Na = -99
	Na = Nw
#	pa = 0.00718
#	sa = 
	pd1 = dfe['pa[0]']
	sd1 = dfe['sa[0]']
	if classes == 2:

		pd2 = 1 - pd1-pa
		sd2 = dfe['Es']
	elif classes == 3:
		pd2 = dfe['pa[1]']
		sd2 = dfe['sa[1]']
		pd3 = 1.0 - pa - pd1 - pd2
		sd3 = dfe['Es']

	uneut = kimura_fixation_prob(0.0, Na)
#	print("uneut: "+ str(uneut))
	ua = kimura_fixation_prob(sa, Na)
	ud1 = kimura_fixation_prob(sd1, Na)
	ud2 = kimura_fixation_prob(sd2, Na)
	if classes == 2:
		usel = pa*ua + pd1*ud1 + pd2*ud2
	if classes ==3:
		ud3 = kimura_fixation_prob(sd3, Na)
		usel = pa*ua + pd1*ud1 + pd2*ud2 + pd3*ud3
	return usel/uneut

def est_dfe_dict(estLine):
	myDat = estLine.split(' ')
	dfeDict = {}
	for i,j in zip(myDat[::2],myDat[1::2]):
		dfeDict[i] = float(j)
	return dfeDict

#N1 100 N2 6 t2 17.8613 N3 80 t3 47.5530 Nw 81.48 b -99.0000 Es -2.695909 f0 0.927200808 f2n 0.019793149 pa[0] 0.200626 sa[0] -0.000261 pa[1] 0.001395 sa[1] 0.169006 L -792213.6530 0.253899999894
def main():
#Parameters for 0-fold sites hard-coded for now but will change once it is running.
## Target dN/dS is 0.23
	parser = argparse.ArgumentParser(description="A li'l program to take DFE-alpha output (est_dfe.out) and calculate the expected fixation probabilities for a specific DFE (discrete class dfe only). There needs to be a class of adaptive mutations in the DFE. ADAPTED FROM PK's CODE FOR THE KIMURA FIXATION PROBABILITY")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The est_dfe.out file that you want to analyse")
	parser.add_argument("-c","--classes", 
		required = True,
		dest = "classes",
		type =int, 
		help = "The number of DFE classes")
	parser.add_argument("-d","--dnds", 
		required = False,
		dest = "dnds",
		type =float, 
		help = "**Optional** The dnds for the class of sites you are looking at. Will be used inthe output file to give you a column of the absolute difference between the predicted and the observed",
		default = 0)
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the resulting output file")
	args = parser.parse_args()
	
	dfe = est_dfe_dict(open(args.input).readlines()[0].strip())
	
	if args.classes==2:
		pA = 'pa[1]'
		sA = 'sa[1]'
	if args.classes==3:
		pA = 'pa[2]'
		sA = 'sa[2]'

	pa = dfe[pA]
	sa = dfe[sA]
	compound = pa*sa
	range_of_pa = [ float(i)/1000000 for i in range(1,100000,1)]
	cor_sa = [ compound/i for i in range_of_pa] # correspoding sA, given pA*sA
	output = [[i, j, i*j, evaluateKaKs(i,j,dfe,args.classes),abs(evaluateKaKs(i,j,dfe,args.classes) - args.dnds)] for i,j in zip(range_of_pa, cor_sa)]

	x = pd.DataFrame(output,columns=['pa','sa','pasa','KaKs','diff'])
	x.to_csv(args.output)

if '__name__':
	main()
