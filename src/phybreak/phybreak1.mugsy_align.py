import os

project_dir = ""
contig_dir = ""
contig_extension = ""
strain_list_file = ""
output_prefix = ""

parameter_file = open("phybreak_parameters.txt","r")
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
		# elif line[0] == "pop_infile_name":
		# 	pop_infile_name = line[1].split(" #")[0]
		# elif line[0] == "ref_iso":
		# 	ref_iso = line[1].split(" #")[0]
		# elif line[0] == "ref_contig":
		# 	ref_contig = line[1].split(" #")[0]
		# elif line[0] == "len_block_threshold":
		# 	len_block_threshold = int(line[1].split(" #")[0])
		# elif line[0] == "gap_prop_thresh":
		# 	gap_prop_thresh = float(line[1].split(" #")[0])
		# elif line[0] == "window_size":
		# 	window_size = int(line[1].split(" #")[0])
parameter_file.close()

output_dir = project_dir+"align/"
slurm_prefix = "#!/bin/bash\n#SBATCH -N 1#SBATCH -n 1\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=40000\n#SBATCH --time=48:00:00\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=username@gmail.com\ncd "+project_dir+"\n"

#################   MAIN   #################
import os
import sys
if os.path.isdir(output_dir) == False:
	os.makedirs(output_dir)

#import os
#lis = [f for f in os.listdir(filedir) if f.endswith(".fa")]
infile = open(project_dir+strain_list_file,"r")
lis = []
for line in infile:
    line = line.strip()
    line = line + contig_extension
    #line = line.replace("-","")
    lis.append(line)
infile.close()

mugsyline = "source /home/username/mugsyenv.sh"+"\n"+"mugsy "
mauveline = "progressiveMauve"
mauveline = mauveline + " --output="+output_dir+output_prefix+".xmfa "

a = 1
for strain in lis:
    mugsyline = mugsyline + contig_dir+strain +" "
    mauveline = mauveline + contig_dir+strain +" "

mugsyline = mugsyline + "--directory "+output_dir+" --prefix "+output_prefix

shfile = open(project_dir + output_prefix +".align.sh","w")
shfile.write(slurm_prefix+ mugsyline+"\n")
#shfile.write(slurm_prefix+mauveline+"\n")
shfile.close()
command = "sbatch "+project_dir + output_prefix +".align.sh"
os.system(command)