import gzip, sys
from site_frequency_spectrum import SFS_from_frequencies as sfs_f

sample_size = sys.argv[2]
output = open(sys.argv[3],"w")
mutation_dict = {}
element_dict ={}
for i in gzip.open(sys.argv[1]):
	x = i.strip("\n").split(" ")
	if x[0].startswith("g"):
		if x[1].startswith("m"):
			continue
		length = int(x[2])-int(x[1])+1
		if x[0] not in element_dict.keys():
			element_dict[x[0]] = length
		else:
			element_dict[x[0]] += length
	try:
		if x[1].startswith("m") and x[4] == "0.5":
			pass
		else: 
			continue
	except IndexError:
		continue
	if x[1] not in mutation_dict.keys():
		mutation_dict[x[1]] = [int(x[7])]
	else:
		mutation_dict[x[1]].append(int(x[7]))


output.write("sample_size\n"+sample_size+"\n")

for i in element_dict.keys():
	output.write(i+"\n"+str(element_dict[i])+"\n")
print mutation_dict
for i in mutation_dict.keys():	
	output.write(i+"\n"+",".join(mutation_dict[i])+"\n")

