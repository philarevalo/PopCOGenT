'''
Produces multiple genome alignment
'''
import os
import sys

project_dir = ""
contig_dir = ""
contig_extension = ""
strain_list_file = ""
output_prefix = ""

with open("phybreak_parameters.txt","r") as paramter_file:
	for line in parameter_file:
		line = line.strip().split(" = ")
		if len(line) > 1:
			if line[0] == "project_dir":
				project_dir = line[1].split(" #")[0]
			elif line[0] == "contig_dir":
				contig_dir = line[1].split(" #")[0]
			elif line[0] == "contig_extension":
				contig_extension = line[1].split(" #")[0]
			elif line[0] == "strain_list_file":
				strain_list_file = line[1].split(" #")[0]
			elif line[0] == "output_prefix":
				output_prefix = line[1].split(" #")[0]

output_dir = project_dir+"align/"
slurm_prefix = "#!/bin/bash\n#SBATCH -N 1#SBATCH -n 1\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=40000\n#SBATCH --time=48:00:00\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=username@gmail.com\ncd "+project_dir+"\n"

#################   MAIN   #################


if os.path.isdir(output_dir) == False:
	os.makedirs(output_dir)

with open(project_dir+strain_list_file) as infile:
	lis = []
	for line in infile:
	    line = line.strip()
	    line = line + contig_extension
	    lis.append(line)

mugsyline = "source /home/username/mugsyenv.sh"+"\n"+"mugsy "

for strain in lis:
    mugsyline = mugsyline + contig_dir+strain +" "

mugsyline = mugsyline + "--directory " + output_dir + " --prefix " + output_prefix

with open(project_dir + output_prefix +".align.sh","w") as shfile:
	shfile.write(slurm_prefix+ mugsyline+"\n")
shfile.close()
command = "sbatch "+project_dir + output_prefix +".align.sh"
os.system(command)