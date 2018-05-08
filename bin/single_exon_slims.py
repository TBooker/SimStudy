import tom,sys,glob
import site_frequency_spectrum as SFS
## little python script to get the SFS for a set of SLIMulations that model a single exon. 75% of sites under selection in the exons
## Can then look at how well PK's program for estimating the strength of deleterious mutations compares to Thanasis'
if len(sys.argv) < 2:
	print "You need to give a directory of slim output files as the input to the program"
	sys.exit()
files = glob.glob(sys.argv[1]+"/*")

neutral= []
selected = []
length = 0
for i in files:	
	x = open(i).read().split("\n")	
	z = tom.slim(x)
	print i
	if z.output == False: continue
	if len(neutral) == 0:
		N = z.sampleN
	print i
	length+=z.length
	for b in z.mutations:
		if float(b[3]) ==0.0:
			neutral.append(int(b[7]))
		elif float(b[3]) < 0:
			selected.append(int(b[7]))
	


nSFS=SFS.SFS_from_frequencies(neutral,length*0.25,N)
sSFS=SFS.SFS_from_frequencies(selected,length*0.75,N)
print "Rename the file, simulated_output.txt as sfs.txt and run in DFE alpha if you like..."
w = open("simulated_output.txt","w")
w.write("1\n"+str(N)+"\n"+" ".join(map(str,sSFS))+"\n"+" ".join(map(str,nSFS)))
